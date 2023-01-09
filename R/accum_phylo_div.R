#' Phylorgenetic Diversity Accumulation of a Community
#' 
#' Diversity and Entropy Accumulation Curves represent the accumulation of 
#' entropy with respect to the sample size.
#' 
#' `accum_phylo_ent()` or `accum_phylo_div()` estimate the phylogenetic 
#' diversity or entropy accumulation curve of a distribution.
#' See [ent_tsallis] for details about the computation of entropy at each level
#' of interpolation and extrapolation.
#' 
#' In accumulation curves, extrapolation if done by estimating the asymptotic 
#' distribution of the community and estimating entropy at different levels 
#' by interpolation. 
#' 
#' Interpolation and extrapolation of integer orders of diversity are from 
#' \insertCite{Chao2014;textual}{divent}.
#' The asymptotic richness is adjusted so that the extrapolated part of the 
#' accumulation joins the observed value at the sample size.
#' 
#' "accumulation" objects can be plotted.
#' They generalize the classical Species Accumulation Curves (SAC) which are 
#' diversity accumulation of order \eqn{q=0}.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities].
#' @param ... Unused.
#'
#' @return A tibble with the site names, the estimators used and the accumulated entropy
#' or diversity at each level of sampling effort.
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' autoplot(accum_phylo_div(paracou_6_abd))
#' 
#' @name accum_phylo_div
NULL


#' @rdname accum_phylo_div
#'
#' @export
accum_phylo_ent <- function(x, ...) {
  UseMethod("accum_phylo_ent")
}


#' @rdname accum_phylo_div
#'
#' @param levels The levels, i.e. the sample sizes of interpolation or 
#' extrapolation: a vector of integer values.
#' 
#' @export
accum_phylo_ent.numeric <- function(
    x,
    tree,
    q = 0,
    normalize = TRUE,
    levels = NULL, 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Set levels if needed
    if (is.null(levels)) {
      sample_size <- max(
        rowSums(
          x[, !colnames(x) %in% non_species_columns]
        )
      )
      levels <- seq_len(sample_size)
    }
  }
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator) 
  
  # TODO
 }


#' @rdname accum_phylo_div
#'
#' @export
accum_phylo_ent.abundances <- function(
    x,
    tree,
    q = 0,
    normalize = TRUE,
    levels = NULL, 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    gamma = FALSE,
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Set levels if needed
    if (is.null(levels)) {
      sample_size <- max(
        rowSums(
          x[, !colnames(x) %in% non_species_columns]
        )
      )
      levels <- seq_len(sample_size)
    }
  }
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator) 
  
  # Calculate abundances along the tree, that are a list of matrices
  if (gamma) {
    the_phylo_abd <- phylo_abd(abundances = metacommunity(x), tree = tree)
  } else {
    the_phylo_abd <- phylo_abd(abundances = x, tree = tree)
  }
  
  # Calculate the entropy
  the_entropy <- vapply(
    levels,
    FUN = function(level) {
      phylo_entropy.phylo_abd(
        phylo_abd = the_phylo_abd,
        tree = tree,
        normalize = normalize,
        # Arguments for ent_tsallis.numeric
        q = q,
        # estimator is not used for interpolation/extrapolation
        estimator = "naive",
        level = level, 
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha, 
        jack_max = jack_max,
        coverage_estimator = coverage_estimator
      )
    },
    FUN.VALUE = 0
  )
  
  
  # TODO

  class(the_accum_phylo_ent) <- c("accumulation", class(the_accum_phylo_ent))
  
  return(the_accum_phylo_ent)
}


#' @rdname accum_phylo_div
#'
#' @export
accum_phylo_div <- function(x, ...) {
  UseMethod("accum_phylo_div")
}


#' @rdname accum_phylo_div
#'
#' @export
accum_phylo_div.numeric <- function(
    x,
    tree,
    q = 0,
    normalize = TRUE,
    levels = NULL, 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Set levels if needed
    if (is.null(levels)) {
      sample_size <- max(
        rowSums(
          x[, !colnames(x) %in% non_species_columns]
        )
      )
      levels <- seq_len(sample_size)
    }
  }
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator) 
  
  # TODO
}


#' @rdname accum_phylo_div
#'
#' @export
accum_phylo_div.abundances <- function(
    x,
    tree,
    q = 0,
    normalize = TRUE,
    levels = NULL, 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    gamma = FALSE,
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Set levels if needed
    if (is.null(levels)) {
      sample_size <- max(
        rowSums(
          x[, !colnames(x) %in% non_species_columns]
        )
      )
      levels <- seq_len(sample_size)
    }
  }
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator) 
  
  # TODO
}
