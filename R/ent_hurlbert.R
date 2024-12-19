#' Hurlbert Entropy of a Community
#'
#' Estimate the Hurlbert entropy \insertCite{Hurlbert1971}{divent} of species from abundance or probability data.
#' Several estimators are available to deal with incomplete sampling.
#'
#' Bias correction requires the number of individuals.
#' See [div_hurlbert] for estimators.
#'
#' Hurlbert's entropy cannot be estimated at a specified level of interpolation or
#' extrapolation, and entropy partitioning is not available.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @seealso [div_hurlbert]
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Entropy of each community
#' ent_hurlbert(paracou_6_abd, k = 2)
#'
#' @name ent_hurlbert
NULL


#' @rdname ent_hurlbert
#'
#' @export
ent_hurlbert <- function(x, k = 2, ...) {
  UseMethod("ent_hurlbert")
}


#' @rdname ent_hurlbert
#'
#' @param estimator An estimator of entropy.
#'
#' @export
ent_hurlbert.numeric <- function(
    x,
    k = 2,
    estimator = c("Hurlbert", "naive"),
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }

  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    the_entropy <- sum(1 - (1 - x)^k)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive",
          order = k,
          entropy = the_entropy
        )
      )
    }
  }

  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)
  if (k > sample_size) {
    stop("The order of diversity cannot be greater than the size of the sample.")
  }
  # Number of observed species
  s_obs <- length(abd)

  # Entropy of a vector of abundances ----
  if (!is_integer_values(abd)) {
    warning("The estimator can't be applied to non-integer values.")
    estimator <- "naive"
  }
  # Naive estimator
  if (estimator == "naive") {
    prob <- abd / sample_size
    the_entropy <- sum(1 - (1 - prob)^k)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive",
          order = k,
          entropy = the_entropy
        )
      )
    }
  }
  # Hurlbert estimator
  if (estimator == "Hurlbert") {
    # Use lchoose and differences to avoid Inf
    lcnk <- lchoose(sample_size, k)
    the_entropy <- s_obs - sum(exp(lchoose(sample_size - abd, k) - lcnk))
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Hurlbert",
          order = k,
          entropy = the_entropy
        )
      )
    }
  }

  warning("estimator was not recognized")
  return(NA)
}



#' @rdname ent_hurlbert
#'
#' @export
ent_hurlbert.species_distribution <- function(
    x,
    k = 2,
    estimator = c("Hurlbert", "naive"),
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }

  # Apply ent_hurlbert.numeric() to each site
  ent_hurlbert_sites <- apply(
    # Eliminate site and weight columns
    x[, !colnames(x) %in% non_species_columns],
    # Apply to each row
    MARGIN = 1,
    FUN = ent_hurlbert.numeric,
    # Arguments
    k = k,
    estimator = estimator,
    as_numeric = as_numeric,
    check_arguments = FALSE
  )
  if (as_numeric) {
    return(ent_hurlbert_sites)
  } else {
    return(
      # Make a tibble with site, estimator and entropy
      tibble::tibble(
        # Restore non-species columns
        x[colnames(x) %in% non_species_columns],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, ent_hurlbert_sites)
      )
    )
  }
}
