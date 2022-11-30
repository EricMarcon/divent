#' Hill number of a Community
#' 
#' Estimate the diversity sensu stricto, i.e. the Hill number of species from abundance or probability data.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#'
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities]
#' @param ... Unused.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#'
#' @return A tibble with the site names, the estimators used and the estimated diversity.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' div_hill(paracou_6_abd, q = 2)
#' div_hill(paracou_6_abd, q = 2, gamma = TRUE)
#' 
#' @name div_hill
NULL


#' @rdname div_hill
#'
#' @export
div_hill <- function(x, q = 1, ...) {
  UseMethod("div_hill")
}


#' @rdname div_hill
#'
#' @param q The order of diversity.
#' @param estimator An estimator of entropy. 
#' @param level The level of interpolation or extrapolation. 
#' It may be a chosen sample size (an integer) or a sample coverage 
#' (a number between 0 and 1). 
#' Richness extrapolation require its asymptotic estimation depending on the 
#' choice of the `estimator`.
#' @param probability_estimator One of the estimators of a probability distribution: 
#' "naive" (the default value), or "Chao2013", "Chao2015", "ChaoShen" to estimate
#' the probabilities of the observed species in the asymptotic distribution.
#' @param unveiling One of the possible unveiling methods to estimate the probabilities 
#' of the unobserved species: "none" (default, no species is added), "uniform" 
#' (all unobserved species have the same probability) or "geometric" (the 
#' unobserved species distribution is geometric).
#' @param richness_estimator An estimator of richness to evaluate the total number of species,
#' see [div_richness].
#' @param jack_alpha The risk level, 5% by default, used to optimize the jackknife order.
#' @param jack_max The highest jackknife order allowed. Default is 10. 
#' @param sample_coverage The sample coverage of `x` calculated elsewhere. 
#' Used to calculate the gamma diversity of meta-communities, see details. 
#' @param as_numeric If `TRUE`, a number is returned rather than a tibble.
#' 
#' @export
div_hill.numeric <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak",
                  "naive"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    sample_coverage = NULL,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  the_entropy <- ent_tsallis.numeric(
    x, 
    q = q, 
    estimator = estimator,
    level = level, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max,
    sample_coverage = sample_coverage,
    as_numeric = TRUE,
    check_arguments = FALSE
  )
  the_diversity <- exp_q(the_entropy, q = q)
  if (as_numeric) {
    return(the_diversity)
  } else {
    return(
      tibble::tibble_row(
        estimator = estimator, 
        order = q,
        diversity = the_diversity
      )
    ) 
  }
}


#' @rdname div_hill
#'
#' @param gamma If `TRUE`, \eqn{\gamma} diversity, i.e. diversity of the metacommunity, is computed.
#' 
#' @export
div_hill.species_distribution <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  if (gamma) {
    # Calculate gamma entropy
    the_entropy <- ent_gamma(
      x = x,
      q = q,
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max
    )
    # Calculate diversity
    the_diversity <- dplyr::mutate(
      the_entropy, 
      diversity = exp_q(.data$entropy, q = q),
      .keep = "unused"
    )
    # return the diversity
    return(the_diversity)
  } else {
    # Apply div_hill.numeric() to each site
    div_hill_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% c("site", "weight"))], 
      # Apply to each row
      MARGIN = 1,
      FUN = div_hill.numeric,
      # Arguments
      q = q,
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max, 
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    return(
      # Make a tibble with site, estimator and richness
      tibble::tibble(
        # Do not assume column site exists
        x[colnames(x) == "site"],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, div_hill_list)
      )
    )
  }
}
