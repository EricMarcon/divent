#' Hill number of a Community
#' 
#' Estimate the diversity sensu stricto, i.e. the \insertCite{Hill1973;textual}{divent} 
#' number of species from abundance or probability data.
#' 
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#' 
#' Estimation techniques are from \insertCite{Chao2003;textual}{divent}, 
#' \insertCite{Grassberger1988;textual}{divent},\insertCite{Holste1998;textual}{divent}, 
#' \insertCite{Bonachela2008;textual}{divent}, \insertCite{Marcon2014a;textual}{divent} 
#' which is actually the max value of "ChaoShen" and "Grassberger", 
#' \insertCite{Zhang2014a;textual}{divent}, \insertCite{Chao2014c;textual}{divent},
#' \insertCite{Chao2015;textual}{divent} and \insertCite{Marcon2015a;textual}{divent}.
#' 
#' The `ChaoJost` estimator \insertCite{Chao2013,Chao2015;textual}{divent} contains 
#' an unbiased part concerning observed species, equal to that of 
#' \insertCite{Zhang2014a;textual}{divent}, and a (biased) estimator of the remaining 
#' bias based on the estimation of the species-accumulation curve. 
#' It is very efficient but slow if the number of individuals is more than a few hundreds.
#' 
#' The unveiled estimators rely on \insertCite{Chao2014c;textual}{divent}, 
#' completed by \insertCite{Marcon2015a;textual}{divent}. 
#' The actual probabilities of observed species are estimated and completed by 
#' a geometric distribution of the probabilities of unobserved species. 
#' The number of unobserved species is estimated by the Chao1 estimator (`UnveilC`), 
#' following \insertCite{Chao2014c;textual}{divent}, or by the iChao1 (`UnveiliC`)
#' or the jackknife (`UnveilJ`).
#' The `UnveilJ` estimator often has a lower bias but a greater variance 
#' \insertCite{Marcon2015a}{divent}.
#' It is a good first choice thanks to the versatility of the jackknife 
#' estimator of richness.
#' 
#' Estimators by \insertCite{Bonachela2008;textual}{divent} and 
#' \insertCite{Holste1998;textual}{divent} are rarely used.
#' 
#' To estimate \eqn{\gamma} diversity, the size of a metacommunity (see 
#' [metacommunity]) is unknown so it has to be set according to a rule which does
#' not ensure that its abundances are integer values.
#' Then, classical bias-correction methods do not apply. 
#' Providing the `sample_coverage` argument allows applying the `ChaoShen` and
#' `Grassberger` estimators to estimate quite well the entropy.
#' 
#' Diversity can be estimated at a specified level of interpolation or 
#' extrapolation, either a chosen sample size or sample coverage 
#' \insertCite{Chao2014}{divent}, rather than its asymptotic value.
#' See [accum_hill] for details.
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
#' div_hill(paracou_6_abd, q = 2)
#' # gamma diversity
#' div_hill(paracou_6_abd, q = 2, gamma = TRUE)
#' 
#' # At 80% coverage
#' div_hill(paracou_6_abd, q = 2, level = 0.8)
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
#' @param estimator an estimator of asymptotic diversity.
#' 
#' @export
div_hill.numeric <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q_threshold = 10,
    sample_coverage = NULL,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator)

  if (q > q_threshold) {
    # Apply the naive estimator of diversity because rounding errors 
    # of exp_q(entropy) are greater than its bias
    prob <- x / sum(x)
    the_diversity <- sum(prob^q)^(1/(1 - q))
    # Possible rounding error again
    if (is.infinite(the_diversity)) {
      # Berger Parker index
      the_diversity <- 1 / max(prob)
    }
    if (as_numeric) {
      return(the_diversity)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive", 
          order = q,
          diversity = the_diversity
        )
      )  
    }
  } else {
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
      coverage_estimator = coverage_estimator,
      sample_coverage = sample_coverage,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    # Calculate diversity
    the_diversity <- dplyr::mutate(
      the_entropy, 
      diversity = exp_q(.data$entropy, q = q),
      .keep = "unused"
    )
    # return the diversity
    if (as_numeric) {
      return(the_diversity$diversity)
    } else {
      return(the_diversity)
    }
  }
}


#' @rdname div_hill
#'
#' @export
div_hill.species_distribution <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q_threshold = 10,
    gamma = FALSE,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator)

  if (q > q_threshold) {
    # Apply the naive estimator of diversity because rounding errors 
    # of exp_q(entropy) are greater than its bias
    div_hill_sites <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns], 
      # Apply to each row
      MARGIN = 1,
      FUN = function(distribution) {
        prob <- distribution / sum(distribution)
        the_diversity <- sum(prob^q)^(1/(1 - q))
        if (is.infinite(the_diversity)) {
          # Berger Parker index if rounding errors
          the_diversity <- 1 / max(prob)
        }
      }
    )
    if (as_numeric) {
      return(div_hill_sites)
    } else {
      return(
        # Make a tibble with site, estimator and diversity
        tibble::tibble(
          # Restore non-species columns
          x[colnames(x) %in% non_species_columns],
          estimator = estimator, 
          order = q,
          diversity = div_hill_sites
        )
      )
    }
  } else {
    # Estimate diversity and transform it into diversity
    the_entropy <- ent_tsallis.species_distribution(
      x, 
      q = q, 
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      gamma = gamma,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
    # Calculate diversity
    if (as_numeric) {
      the_diversity = exp_q(the_entropy, q = q)
    } else {
      the_diversity <- dplyr::mutate(
        the_entropy, 
        diversity = exp_q(.data$entropy, q = q),
        .keep = "unused"
      )
    }
  }
  return(the_diversity)
}
