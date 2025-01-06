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
#' @param x An object, that may be a named numeric vector (names are species names)
#' containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated
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

  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
    # Prepare the tree
    tree <- as_phylo_divent(tree)
    # Check species names
    species_names <- names(x)
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      cli::cli_abort("Some species are missing in the tree.")
    }
  }

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
    as_numeric = as_numeric,
    check_arguments = FALSE)

  # Return
  return(the_entropy)
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
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
    # Prepare the tree
    tree <- as_phylo_divent(tree)
    # Check species names
    col_names <- colnames(x)
    species_names <- col_names[!col_names %in% non_species_columns]
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      cli::cli_abort("Some species are missing in the tree.")
    }
  }

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
  if (as_numeric) {
    return(the_entropy)
  } else {
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
}


#' Phylogenetic entropies
#'
#' Calculate entropies of a list of phylogenetic abundances (obtained by
#' [phylo_abd]).
#' Each item of the list corresponds to a phylogenetic group, i.e. an interval
#' of the tree (where the species do not change).
#'
#' @param phylo_abd A list of matrices of abundance (caution: rows are species,
#' columns are communities).
#'
#' @returns A vector. Each item is the entropy of a community.
#'
#' @noRd
#'
phylo_entropy.phylo_abd <- function(
  q,
  phylo_abd,
  tree,
  normalize,
  # Other arguments for ent_tsallis.numeric
  estimator,
  level,
  probability_estimator,
  unveiling,
  richness_estimator,
  jack_alpha,
  jack_max,
  coverage_estimator,
  gamma) {

  if (gamma) {
    # Calculate gamma entropy of each group.
    # simplify2array() makes a vector with the list of numbers.
    phylo_entropies <- simplify2array(
      lapply(
        # Calculate entropy in each item of the list, i.e. group.
        # Obtain a list.
        phylo_abd,
        FUN = function(group) {
          ent_tsallis.species_distribution(
            as_abundances.numeric(t(group)),
            q = q,
            estimator = estimator,
            level = level,
            probability_estimator = probability_estimator,
            unveiling = unveiling,
            richness_estimator = richness_estimator,
            jack_alpha  = jack_alpha,
            jack_max = jack_max,
            coverage_estimator = coverage_estimator,
            gamma = TRUE,
            # Obtain a vector.
            as_numeric = TRUE,
            check_arguments = FALSE
          )
        }
      )
    )

  } else {
    # Calculate entropy of each community in each group.
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
            # Arguments
            q = q,
            estimator = estimator,
            level = level,
            probability_estimator = probability_estimator,
            unveiling = unveiling,
            richness_estimator = richness_estimator,
            jack_alpha  = jack_alpha,
            jack_max = jack_max,
            coverage_estimator = coverage_estimator,
            # Obtain a vector.
            as_numeric = TRUE,
            check_arguments = FALSE
          )
        }
      )
    )
  }
  # Should be a matrix, but simplify2array() makes a vector instead of a 1-col
  # matrix and gamma entropy is a vector. Force a matrix.
  if (is.vector(phylo_entropies)) {
    phylo_entropies <- t(phylo_entropies)
  }

  # Calculate the weighted mean of entropy and normalize
  the_entropy <- as.numeric(tree$intervals %*% t(phylo_entropies))
  if (normalize) the_entropy <- the_entropy / sum(tree$intervals)

  return(the_entropy)
}
