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
#' See [accum_tsallis] for details.
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

  # Make a species_distribution
  species_distribution <- as_species_distribution(x)
  
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
    coverage_estimator = coverage_estimator,
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
  
  # Calculate abundances along the tree, that are a list of matrices
  the_phylo_abd <- phylo_abd(abundances = x, tree = tree)

  # Calculate the entropy
  the_entropy <- phylo_entropy.phylo_abd(
    phylo_abd = the_phylo_abd,
    tree = tree,
    normalize = normalize,
    # Arguments for ent_tsallis.numeric
    q = q,
    estimator = estimator,
    level = level, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    gamma = gamma
  )
  
  # Return
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
