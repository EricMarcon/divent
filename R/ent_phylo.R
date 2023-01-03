#' Phylogenetic Entropy of a Community
#' 
#' Estimate the entropy of species from abundance or probability data
#' and a phylogenetic tree.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#' See [div_hill] for estimators.
#' 
#' Entropy can be estimated at a specified level of interpolation or 
#' extrapolation, either a chosen sample size or sample coverage 
#' \insertCite{Chao2014}{divent}, rather than its asymptotic value.
#' See [ent_accum] for details.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances 
#' or probabilities, or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @return A tibble with the site names, the estimators used and the estimated
#' entropy.
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Entropy of each community
#' ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2)
#' # Gamma entropy
#' ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2, gamma = TRUE)
#' 
#' # At 80% coverage
#' ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2, level = 0.8)
#' 
#' 
#' @name ent_phylo
NULL


#' @rdname ent_phylo
#'
#' @export
ent_phylo <- function(x, tree, q = 1, ...) {
  UseMethod("ent_phylo")
}


#' @rdname ent_phylo
#'
#' @param estimator An estimator of entropy. 
#' 
#' @export
ent_phylo.numeric <- function(
    x,
    tree,
    q = 1,
    normalize = TRUE,
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  tree <- as_phylo_divent(tree)

  # Make a species_distribution
  species_distribution <- as_species_distribution(as.vector(x))
  
  # Entropy
  the_entropy <- ent_phylo.species_distribution(
    species_distribution,
    tree = tree,
    q = q,
    normalize = normalize,
    estimator = estimator,
    level = level, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max,
    check_arguments = FALSE)
 
  # Return
  if(as_numeric) {
    return(the_entropy$entropy)
  } else {
    return(the_entropy)
  }
}

#' @rdname ent_phylo
#'
#' @export
ent_phylo.species_distribution <- function(
    x,
    tree,
    q = 1,
    normalize = TRUE,
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  tree <- as_phylo_divent(tree)
  
  # Species names
  col_names <- colnames(x)
  species_names <- col_names[!(col_names %in% non_species_columns)]
  
  # Calculate abundances along the tree, that are a list of matrices
  phylo_abd <- sapply(
    # Each phylogenetic group yields an item of the list
    colnames(tree$phylo_groups), 
    function(group) {
      # Create a matrix with the abundances of groups in each community
      apply(
        x[, !(colnames(x) %in% non_species_columns)], 
        # Each community yields a column of the matrix
        MARGIN = 1, 
        FUN = function(abd) {
          tapply(
            as.numeric(abd), 
            # Each group yields a row of the matrix
            INDEX = tree$phylo_groups[names(abd), group], 
            FUN = sum
          )
        }
      )
    }, 
    simplify = FALSE
  )
  
  # Calculate entropy of each community in each group.
  if (gamma) {
    phylo_entropies <- simplify2array(
      lapply(
        # Calculate entropy in each item of the list, i.e. group.
        # Obtain a list.
        phylo_abd,
        FUN = function(group) {
          ent_gamma.matrix(
            abd = group,
            weights = x$weight,
            q = q,
            estimator = estimator,
            level = level, 
            probability_estimator = probability_estimator,
            unveiling = unveiling,
            richness_estimator = richness_estimator,
            jack_alpha  = jack_alpha, 
            jack_max = jack_max,
            coverage_estimator = coverage_estimator,
          )
        }
      )
    )
  } else {
    # simplify2array() makes a matrix with the list of vectors.
    phylo_entropies <- simplify2array(
      lapply(
        # Calculate entropy in each item of the list, i.e. group.
        # Obtain a list.
        phylo_abd,
        FUN = function(group) {
          apply(
            group,
            # Calculate entropy of each column of the matrix, i.e. community.
            MARGIN = 2,
            FUN = ent_tsallis.numeric,
            q = q,
            estimator = estimator,
            level = level, 
            probability_estimator = probability_estimator,
            unveiling = unveiling,
            richness_estimator = richness_estimator,
            jack_alpha  = jack_alpha, 
            jack_max = jack_max,
            # Obtain a vector.
            as_numeric = TRUE,
            check_arguments = FALSE
          )
        }
      )
    )
  }
  # Should be a matrix, but simplify2array() makes a vector instead of a 1-col 
  # matrix. Force a matrix.
  if(is.vector(phylo_entropies)) phylo_entropies <- t(phylo_entropies)
  
  # Calculate the weighted mean of entropy and normalize
  the_entropy <- as.numeric(tree$intervals %*% t(phylo_entropies))
  if (normalize) the_entropy <- the_entropy / sum(tree$intervals)
  
  if (gamma) {
    return(
      # Make a tibble with site, estimator and entropy
      tibble::tibble_row(
        site = "Metacommunity",
        # estimator and order
        estimator = estimator,
        q = q,
        # Entropy
        entropy = the_entropy
      )
    )
  } else {
    return(
      # Make a tibble with site, estimator and entropy
      tibble::tibble(
        # Restore non-species columns
        x[colnames(x) %in% non_species_columns],
        # estimator and order
        estimator = estimator,
        q = q,
        # Entropy
        entropy = the_entropy
      )
    )
  }
}

#' Gamma entropy of a metacommunity
#' 
#' `abd` is assumed to be a matrix of abundances, lines are communities.
#' `weights` are necessary for gamma diversity.
#' 
#' @param abd A matrix containing abundances or probabilities.
#' 
#' @return A number equal to gamma entropy.
#' @noRd
ent_gamma.matrix <- function(
    abd,
    weights,
    q,
    estimator,
    level,
    probability_estimator,
    unveiling,
    richness_estimator,
    jack_alpha,
    jack_max,
    coverage_estimator) {
  
  # Build the species distribution
  species_distribution <- species_distribution(
    t(abd),
    weights = weights,
    check_arguments = FALSE
  )
  
  # Call ent_gamma.species_distribution
  return(
    ent_gamma.species_distribution(
      species_distribution,
      q = q,
      estimator = estimator,
      level = level,
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha = jack_alpha,
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      as_numeric = TRUE
    )
  )
}
