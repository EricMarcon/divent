#' Phylogenetic Diversity Accumulation of a Community
#'
#' Diversity and Entropy Accumulation Curves represent the accumulation of
#' entropy with respect to the sample size.
#'
#' `accum_ent_phylo()` or `accum_div_phylo()` estimate the phylogenetic
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
#' @returns A tibble with the site names, the estimators used and the accumulated entropy
#' or diversity at each level of sampling effort.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Richness accumulation up to the sample size.
#' # 100 simulations only to save time.
#' autoplot(
#'   accum_div_phylo(mock_3sp_abd, tree = mock_3sp_tree, n_simulations = 100)
#' )
#'
#' @name accum_div_phylo
NULL


#' @rdname accum_div_phylo
#'
#' @export
accum_ent_phylo <- function(x, ...) {
  UseMethod("accum_ent_phylo")
}


#' @rdname accum_div_phylo
#'
#' @param levels The levels, i.e. the sample sizes of interpolation or
#' extrapolation: a vector of integer values.
#'
#' @export
accum_ent_phylo.numeric <- function(
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

  # Check arguments
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
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
    # Set levels if needed
    if (is.null(levels)) {
      sample_size <- sum(x)
      levels <- seq_len(sample_size)
    }
  }

  # Make a species_distribution
  the_species_distribution <- as_species_distribution(x)

  # Entropy accumulation
  the_entropy <- accum_ent_phylo.abundances(
    x = the_species_distribution,
    tree = tree,
    q = q,
    normalize = normalize,
    levels = levels,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    gamma = FALSE,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )

  # Return
  return(the_entropy)
}


#' @rdname accum_div_phylo
#'
#' @export
accum_ent_phylo.abundances <- function(
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

  # Check arguments
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    # Prepare the tree
    tree <- as_phylo_divent(tree)
    # Check species names
    if (length(setdiff(species_names, rownames(tree$phylo_groups))) != 0) {
      stop("Some species are missing in the tree.")
    }
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

  # Species names
  col_names <- colnames(x)
  species_names <- col_names[!col_names %in% non_species_columns]

  # Calculate abundances along the tree, that are a list of matrices
  if (gamma) {
    the_phylo_abd <- phylo_abd(abundances = metacommunity(x), tree = tree)
  } else {
    the_phylo_abd <- phylo_abd(abundances = x, tree = tree)
  }

  # Prepare arrays to store entropy (3 dimensions: x, y, z)
  # and simulated entropies (4 dimensions : x, y, z, t)
  # x are tree intervals, y are communities, z are levels,
  # t are simulations.
  # Add an array to store simulation envelopes, where
  # t are quantiles of simulations, inf and sup.
  if (gamma) {
    n_communities <- 1
  } else {
    n_communities <- nrow(x)
  }
  ent_phylo_abd <- array(
    dim = c(length(tree$intervals), nrow(x), length(levels))
  )
  if (n_simulations > 0) {
    ent_phylo_sim <- array(
      dim = c(length(tree$intervals), nrow(x), length(levels), n_simulations)
    )
    ent_phylo_envelope <- array(
      dim = c(length(tree$intervals), nrow(x), length(levels), 2)
    )
  }
  # Prepare the progress bar
  if (show_progress & interactive()) {
    cli::cli_progress_bar(
      "Computing entropy",
      total = (length(the_phylo_abd) * n_communities) * (1 + n_simulations))
  }

  # Calculate entropy along the tree
  for (x_interval in seq_along(the_phylo_abd)) {
    if (n_simulations > 0) {
      # Simulate communities
      comm_sim.list <- apply(
        # Produce a list of abundances, each of them contains n_simulations of
        # a community
        the_phylo_abd[[x_interval]],
        MARGIN = 2,
        FUN = function(abd) {
          rcommunity(
            n = n_simulations,
            abd = abd,
            check_arguments = FALSE
          )
        }
      )
      # Prepare an array to store simulated abd. Rows are species, columns are
      # communities (same structure as groups of the_phylo_abd), z are simulations
      # Max number of species in simulated communities
      sp_sim <- max(
        vapply(
          comm_sim.list,
          FUN = dim,
          FUN.VALUE = c(0L, 0L)
        )[2, ]
      )
      comm_sim <- array(
        data = 0,
        dim = c(sp_sim, length(comm_sim.list), n_simulations)
      )
      # Move the simulations from the list to the array
      for (simulation in seq_len(n_simulations)) {
        for (community in seq_along(comm_sim.list)) {
          # Number of species in the simulation
          # (= number of columns - 2 for site and weight)
          sp_sim <- dim(comm_sim.list[[community]])[2] - 2
          # Pick a simulation. Store simulated species, let extra cols = 0
          comm_sim[1:sp_sim, community, simulation] <- as.numeric(
            # Corresponding item in the list, remove site and weight
            comm_sim.list[[community]][simulation, -(1:2)]
          )
        }
      }
    }

    for (y_community in seq_len(n_communities)) {
      # Calculate the profile of each community
      # Actual data
      ent_phylo_abd[x_interval, y_community, ] <- accum_tsallis.numeric(
        x = the_phylo_abd[[x_interval]][, y_community],
        q = q,
        levels = levels,
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha = jack_alpha,
        jack_max = jack_max,
        coverage_estimator = coverage_estimator,
        n_simulations = 0,
        alpha = alpha,
        show_progress = FALSE,
        check_arguments = FALSE
      )$entropy

      for (t_simulation in seq_len(n_simulations)) {
        # Entropy of simulated communities
        ent_phylo_sim[x_interval, y_community, , t_simulation] <- accum_tsallis.numeric(
          x = comm_sim[, y_community, t_simulation],
          q = q,
          levels = levels,
          probability_estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator,
          jack_alpha = jack_alpha,
          jack_max = jack_max,
          coverage_estimator = coverage_estimator,
          n_simulations = 0,
          alpha = alpha,
          show_progress = FALSE,
          check_arguments = FALSE
        )$entropy
        # Progress bar
        if (show_progress & interactive()) cli::cli_progress_update()
      }
      if (n_simulations > 0) {
        for (y_community in seq_len(n_communities)) {
          for (z_level in seq_along(levels)) {
            # Quantiles, recentered
            ent_phylo_envelope[x_interval, y_community, z_level, ] <- stats::quantile(
              ent_phylo_sim[x_interval, y_community, z_level, ],
              probs = c(alpha / 2, 1 - alpha / 2),
              na.rm = TRUE
            ) - mean(ent_phylo_sim[x_interval, y_community, z_level, ]) +
              ent_phylo_abd[x_interval, y_community, z_level]
          }
        }
      }

      # Progress bar
      if (show_progress & interactive()) cli::cli_progress_update()
    }
  }
  if (show_progress & interactive()) cli::cli_progress_done()

  # Average entropy
  # Actual data
  ent_community <- apply(
    ent_phylo_abd,
    MARGIN = 2:3,
    FUN = stats::weighted.mean,
    # Arguments
    w = tree$intervals
  )
  # Simulations
  if (n_simulations > 0) {
    ent_quantiles <- apply(
      ent_phylo_envelope,
      MARGIN = 2:4,
      FUN = stats::weighted.mean,
      # Arguments
      w = tree$intervals
    )
  }

  # Format the result
  the_profile_phylo <- ent.tibble(
    ent.matrix = ent_community,
    x = x,
    levels = levels
  )
  # Add the estimator
  sample_sizes <- rowSums(x[, species_names])
  names(sample_sizes) = x$site
  the_profile_phylo <- dplyr::mutate(
    the_profile_phylo,
    estimator = dplyr::case_when(
      .data$level == sample_sizes[.data$site] ~ "Sample",
      .data$level < sample_sizes[.data$site] ~ "Interpolation",
      TRUE ~ "Extrapolation"
    ),
    .before = "entropy"
  )
  # Add simulation columns
  if (n_simulations > 0) {
    ent_inf <- ent.tibble(
      ent.matrix = ent_quantiles[, , 1],
      x = x,
      levels = levels
    )
    ent_sup <- ent.tibble(
      ent.matrix = ent_quantiles[, , 2],
      x = x,
      levels = levels
    )
    the_profile_phylo <- dplyr::bind_cols(
      the_profile_phylo,
      inf = ent_inf$entropy,
      sup = ent_sup$entropy
    )
  }

  class(the_profile_phylo) <- c("accumulation", class(the_profile_phylo))
  return(the_profile_phylo)
}


#' @rdname accum_div_phylo
#'
#' @export
accum_div_phylo <- function(x, ...) {
  UseMethod("accum_div_phylo")
}


#' @rdname accum_div_phylo
#'
#' @export
accum_div_phylo.numeric <- function(
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

  # Check arguments
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
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
    # Set levels if needed
    if (is.null(levels)) {
      sample_size <- sum(x)
      levels <- seq_len(sample_size)
    }
  }

  # Make a the_species_distribution
  the_species_distribution <- as_species_distribution(x)

  # Diversity accumulation
  the_diversity <- accum_div_phylo.abundances(
    x = the_species_distribution,
    tree = tree,
    q = q,
    normalize = normalize,
    levels = levels,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    gamma = FALSE,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )

  # Return
  return(the_diversity)
}


#' @rdname accum_div_phylo
#'
#' @export
accum_div_phylo.abundances <- function(
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

  # Check arguments
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
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

  the_entropy <- accum_ent_phylo.abundances(
    x,
    tree = tree,
    q = q,
    normalize = normalize,
    levels = levels,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    gamma = gamma,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )

  # Calculate diversity
  the_diversity <- dplyr::mutate(
    the_entropy,
    diversity = exp_q(.data$entropy, q = q),
    .keep = "unused"
  )
  if (n_simulations > 0) {
    the_diversity <- dplyr::mutate(
      the_diversity,
      inf = exp_q(.data$inf, q = q),
      sup = exp_q(.data$sup, q = q)
    )
  }
  return(the_diversity)
}


#' Make a long tibble of entropy with a matrix of entropy
#'
#' Utility for [accum_ent_phylo.abundances]
#'
#' @param ent.matrix The matrix of entropies.
#' Rows are communities, columns are orders of entropy.
#' @param x The species distribution.
#' @param levels The levels of interpolation and extrapolation.
#'
#' @returns A tibble. Columns are "site", "level" and "entropy".
#' @noRd
#'
ent.tibble <- function(ent.matrix, x, levels) {

  if (!is.matrix(ent.matrix)) {
    # ent.matrix may be a numeric vector (single community / min and max)
    ent.matrix <- t(as.matrix(ent.matrix))
  }
  # Make a tibble with site names and entropies.
  # Columns are levels of inter/extrapolation
  the_ent.tibble <- tibble::tibble(
    site = x$site,
    data.frame(ent.matrix)
  )
  colnames(the_ent.tibble)[-1] <- as.character(levels)
  # Make a long tibble with an "level" column
  the_ent.tibble <- tidyr::pivot_longer(
    the_ent.tibble,
    cols = !.data$site,
    names_to = "level",
    names_transform = list(level = as.numeric),
    values_to = "entropy"
  )
  # Return
  return(the_ent.tibble)
}
