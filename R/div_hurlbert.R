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
#' @return A tibble with the site names, the estimators used and the estimated diversity.
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

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  estimator <- match.arg(estimator) 

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
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  estimator <- match.arg(estimator) 

  the_entropy <- ent_hurlbert.species_distribution(
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

  return(the_diversity)
}
