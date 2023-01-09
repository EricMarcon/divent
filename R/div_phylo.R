#' Phylogenetic Diversity of a Community
#' 
#' Estimate the diversity of species from abundance or probability data
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
#' diversity
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' div_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2)
#' 
#' # At 80% coverage
#' div_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2, level = 0.8)
#' 
#' # Gamma entropy
#' div_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2, gamma = TRUE)
#' 
#' @name div_phylo
NULL


#' @rdname div_phylo
#'
#' @export
div_phylo <- function(x, tree, q = 1, ...) {
  UseMethod("div_phylo")
}


#' @rdname div_phylo
#'
#' @param estimator An estimator of asymptotic diversity.
#' 
#' @export
div_phylo.numeric <- function(
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
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Prepare the tree
    tree <- as_phylo_divent(tree)
    # Check species names
    col_names <- colnames(x)
    species_names <- col_names[!col_names %in% non_species_columns]
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      stop("Some species are missing in the tree.")    
    }
  }
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator)

  the_entropy <- ent_phylo.numeric(
    x, 
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
    coverage_estimator = coverage_estimator,
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
  return(the_diversity)
}


#' @rdname div_phylo
#'
#' @export
div_phylo.species_distribution <- function(
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
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Prepare the tree
    tree <- as_phylo_divent(tree)
    # Check species names
    col_names <- colnames(x)
    species_names <- col_names[!col_names %in% non_species_columns]
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      stop("Some species are missing in the tree.")    
    }
  }
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator)

  the_entropy <- ent_phylo.species_distribution(
    x,
    tree = tree,
    normalize = TRUE,
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
    check_arguments = FALSE
  )
  # Calculate diversity
  the_diversity <- dplyr::mutate(
    the_entropy, 
    diversity = exp_q(.data$entropy, q = q),
    .keep = "unused"
  )

  return(the_diversity)
}
