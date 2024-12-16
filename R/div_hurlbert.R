#' Hurlbert Diversity of a Community
#'
#' Estimate the diversity sensu stricto, i.e. the effective number of species
#' number of species \insertCite{Dauby2012;textual}{divent}
#' from abundance or probability data.
#'
#' Several estimators are available to deal with incomplete sampling.
#'
#' Bias correction requires the number of individuals.
#'
#' Estimation techniques are from \insertCite{Hurlbert1971;textual}{divent}.
#'
#' Hurlbert's diversity cannot be estimated at a specified level of interpolation or
#' extrapolation, and diversity partioning is not available.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities].
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Diversity of each community
#' div_hurlbert(paracou_6_abd, k = 2)
#'
#' @name div_hurlbert
NULL


#' @rdname div_hurlbert
#'
#' @export
div_hurlbert <- function(x, k = 1, ...) {
  UseMethod("div_hurlbert")
}


#' @rdname div_hurlbert
#'
#' @param estimator An estimator of asymptotic diversity.
#'
#' @export
div_hurlbert.numeric <- function(
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

  the_entropy <- ent_hurlbert.numeric(
    x,
    k = k,
    estimator = estimator,
    check_arguments = FALSE
  )
  # Calculate diversity
  the_diversity <- dplyr::mutate(
    the_entropy,
    diversity = hurlbert_ent2div(.data$entropy, k = k),
    .keep = "unused"
  )

  # return the diversity
  if (as_numeric) {
    return(the_diversity$diversity)
  } else {
    return(the_diversity)
  }
}


#' @rdname div_hurlbert
#'
#' @export
div_hurlbert.species_distribution <- function(
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

  the_entropy <- ent_hurlbert.species_distribution(
    x,
    k = k,
    estimator = estimator,
    as_numeric = as_numeric,
    check_arguments = FALSE
  )
  # Calculate diversity
  if (as_numeric) {
    the_diversity = hurlbert_ent2div(the_entropy, k = k)
  } else {
    the_diversity <- dplyr::mutate(
      the_entropy,
      diversity = hurlbert_ent2div(.data$entropy, k = k),
      .keep = "unused"
    )
  }

  return(the_diversity)
}


#' Compute Hurlbert's diversity from its entropy
#'
#' Find the effective number of species numerically
#'
#' @param hurlbert_entropy The entropy.
#' @param k The order of entropy.
#'
#' @returns Hurlbert's effective number of species.
#' @noRd
#'
hurlbert_ent2div <- function(hurlbert_entropy, k) {
  # Relation between diversity and entropy
  # (D for diversity, S for entropy, k is the parameter)
  f <- function(D, S, k) {D * (1 - (1 - 1 / D)^k) - S}
  # Minimize it
  return(
    vapply(
      hurlbert_entropy,
      FUN = function(S) {
        stats::uniroot(
          f = f,
          interval = c(1, 1E+7),
          S = S,
          k = k
        )$root
      },
      FUN.VALUE = 0
    )
  )
}
