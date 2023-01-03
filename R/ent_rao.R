#' Rao's Quadratic Entropy of a Community
#' 
#' Estimate the quadratic entropy \insertCite{Rao1982}{divent} of species from 
#' abundance or probability data.
#' An estimator \insertCite{Lande1996}{divent} is available to deal with 
#' incomplete sampling.
#' 
#' Rao's entropy is phylogenetic or similarity-based entropy of order 2.
#' [ent_phylo] and [ent_similarity] with argument `q = 2` provide more estimators
#' and allow estimating entropy at a chosen level.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @return A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Entropy of each community
#' ent_rao(paracou_6_abd, tree = paracou_6_taxo)
#' # Similar to (but estimators are not the same) 
#' ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2)
#' 
#' # Entropy of each community
#' ent_rao(paracou_6_abd, distances = paracou_6_fundist)
#' # Similar to (but estimators are not the same) 
#' ent_similarity(
#'   paracou_6_abd, 
#'   similarities = fun_similarity(paracou_6_fundist), 
#'   q = 2
#' )
#' 
#' # gamma entropy
#' ent_rao(paracou_6_abd, tree = paracou_6_taxo, gamma = TRUE)
#' 
#' @name ent_rao
NULL


#' @rdname ent_rao
#'
#' @export
ent_rao <- function(x, ...) {
  UseMethod("ent_rao")
}


#' @rdname ent_rao
#'
#' @param estimator An estimator of entropy.
#' 
#' @export
ent_rao.numeric <- function(
    x, 
    distances = NULL,
    tree = NULL,
    normalize = TRUE,
    estimator = c("Lande", "naive"),
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  if (!xor(is.null(distances), is.null(tree))) {
    stop("Either 'distance' or 'tree' must be provided.")
  }
  
  # Calculate distances from tree
  if (is.null(distances)) {
    tree <- as_phylo_divent(tree)
    distances <- as.matrix(stats::cophenetic(tree$hclust))
  }
  # Check that dimensions correspond
  if (length(x) != ncol(distances)) {
    stop("The length of 'x' must be equal to the dimension of the distances.")
  }
  # Normalize
  if (normalize) {
    distances <- distances / max(distances)
  }
  
  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    if (!is.null(level)) stop("Entropy can't be estimated at a level without the abundance of species.")
    # Probabilities sum to 1, allowing rounding error
    the_entropy <- mean(distances %*% x)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive", 
          order = 2,
          entropy = the_entropy
        )
      )  
    }
  }
  
  # Eliminate 0
  abd <- as.numeric(x[x > 0])
  distances <- distances[x > 0, x > 0]
  # Sample size
  sample_size <- sum(abd)
  
  ## Exit if x contains no or a single species ----
  if (length(abd) < 2) {
    if (length(abd) == 0) {
      if (as_numeric) {
        return(NA)
      } else {
        return(
          tibble::tibble_row(
            estimator = "No Species", 
            order = 2,
            entropy = NA
          )
        )  
      }
    } else {
      if (as_numeric) {
        return(0)
      } else {
        return(
          tibble::tibble_row(
            estimator = "Single Species", 
            order = 2,
            entropy = 0
          )
        )  
      }
    }
  } else {
    # Probabilities instead of abundances
    if (sample_size < 2) {
      warning("Entropy estimators can't apply to probability data. Estimator forced to 'naive'")
      estimator <- "naive"
    }
  }
  
  # Probabilities
  prob <- abd / sum(abd)
  
  if (estimator == "naive") {
    # Naive estimator ----
    the_entropy <- mean(distances %*% prob)
  } else if (estimator == "Lande") {
    # Lande's estimator ----
    the_entropy <- mean(distances %*% prob) * sample_size / (sample_size - 1)
  }
  
  # Return
  if (as_numeric) {
    return(the_entropy)
  } else {
    return(
      tibble::tibble_row(
        estimator = estimator, 
        order = 2,
        entropy = the_entropy
      )
    )  
  }     
}


#' @rdname ent_rao
#'
#' @export
ent_rao.species_distribution <- function(
    x,
    distances = NULL,
    tree = NULL,
    normalize = TRUE,
    estimator = c("Lande", "naive"),
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  if (!xor(is.null(distances), is.null(tree))) {
    stop("Either 'distance' or 'tree' must be provided.")
  }
  
  if (gamma) {
    # Build the metacommunity
    abd <- metacommunity(
      x, 
      as_numeric = TRUE, 
      check_arguments = FALSE
    )
    the_entropy <- ent_rao.numeric(
      abd,
      distances = distances,
      tree = tree,
      normalize = normalize,
      estimator = estimator,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    # Add the site column
    the_entropy <- dplyr::bind_cols(
      site = "Metacommunity",
      the_entropy
    )
    return(the_entropy)
  } else {
    # Apply ent_rao.numeric() to each site
    ent_rao_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% non_species_columns)], 
      # Apply to each row
      MARGIN = 1,
      FUN = ent_rao.numeric,
      # Arguments
      distances = distances,
      tree = tree,
      normalize = normalize,
      estimator = estimator,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    return(
      # Make a tibble with site, estimator and richness
      tibble::tibble(
        # Restore non-species columns
        x[colnames(x) %in% non_species_columns],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, ent_rao_list)
      )
    )
  }
}
