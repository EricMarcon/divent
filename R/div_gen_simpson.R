#' Generalized Simpson's Diversity
#'
#' Estimate the diversity sensu stricto, i.e. the effective number of species
#' \insertCite{Grabchak2016}{divent} from abundance or probability data.
#'
#' Bias correction requires the number of individuals.
#'
#' Estimation techniques are from \insertCite{Zhang2014;textual}{divent}.
#' It is limited to orders \eqn{k} less than or equal to the number of individuals
#' in the community.
#'
#' Generalized Simpson's diversity cannot be estimated at a specified level
#' of interpolation or extrapolation, and diversity partitioning is not available.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities].
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity.
#'
#' @seealso [ent_gen_simpson]
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Diversity of each community
#' div_gen_simpson(paracou_6_abd, k = 50)
#'
#' @name div_gen_simpson
NULL


#' @rdname div_gen_simpson
#'
#' @export
div_gen_simpson <- function(x, k = 1, ...) {
  UseMethod("div_gen_simpson")
}


#' @rdname div_gen_simpson
#'
#' @param estimator An estimator of asymptotic diversity.
#'
#' @export
div_gen_simpson.numeric <- function(
    x,
    k = 1,
    estimator = c("Zhang", "naive"),
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }

  the_entropy <- ent_gen_simpson.numeric(
    x,
    k = k,
    estimator = estimator,
    as_numeric = FALSE,
    check_arguments = FALSE
  )
  # Calculate diversity
  the_diversity <- dplyr::mutate(
    the_entropy,
    diversity = 1 / (1 - .data$entropy)^(1 / k),
    .keep = "unused"
  )

  # return the diversity
  if (as_numeric) {
    return(the_diversity$diversity)
  } else {
    return(the_diversity)
  }
}


#' @rdname div_gen_simpson
#'
#' @export
div_gen_simpson.species_distribution <- function(
    x,
    k = 1,
    estimator = c("Zhang", "naive"),
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }

  the_entropy <- ent_gen_simpson.species_distribution(
    x,
    k = k,
    estimator = estimator,
    as_numeric = as_numeric,
    check_arguments = FALSE
  )
  # Calculate diversity
  if (as_numeric) {
    the_diversity <- 1 / (1 - the_entropy)^(1 / k)
  } else {
    the_diversity <- dplyr::mutate(
      the_entropy,
      diversity = 1 / (1 - .data$entropy)^(1 / k),
      .keep = "unused"
    )
  }

  return(the_diversity)
}
