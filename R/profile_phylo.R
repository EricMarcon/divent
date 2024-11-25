#' Phylogenetic Diversity Profile of a Community
#'
#' Calculate the diversity profile of a community, i.e. its phylogenetic diversity
#' against its order.
#'
#' A bootstrap confidence interval can be produced by simulating communities
#' (their number is `n_simulations`) with [rcommunity] and calculating their profiles.
#' Simulating communities implies a downward bias in the estimation:
#' rare species of the actual community may have abundance zero in simulated communities.
#' Simulated diversity values are recentered so that their mean is that of the actual community.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @examples
#' profile_phylo(paracou_6_abd, tree = paracou_6_taxo)
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity at each order.
#' This is an object of class "profile" that can be plotted.
#'
#' @references
#' \insertAllCited{}
#'
#' @name profile_phylo
NULL


#' @rdname profile_phylo
#'
#' @export
profile_phylo <- function(
    x,
    tree,
    orders = seq(from = 0, to = 2, by = 0.1),
    ...) {
  UseMethod("profile_phylo")
}


#' @rdname profile_phylo
#'
#' @param orders The orders of diversity used to build the profile.
#' @param estimator An estimator of entropy.
#' @param n_simulations The number of simulations used to estimate the confidence envelope of the profile.
#' @param alpha The risk level, 5% by default, of the confidence envelope of the profile.
#'
#' @export
profile_phylo.numeric <- function(
    x,
    tree,
    orders = seq(from = 0, to = 2, by = 0.1),
    normalize = TRUE,
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger",
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak",
                  "naive"),
    level = NULL,
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05,
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    sample_coverage = NULL,
    as_numeric = FALSE,
    n_simulations = 0,
    alpha = 0.05,
    bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
    show_progress = TRUE,
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
  bootstrap <- match.arg(bootstrap)
  if (as_numeric && n_simulations > 0) {
    stop ("No simulations are allowed if a numeric vector is expected ('as_numeric = TRUE').")
  }

  # Call the .species_distribution method
  the_profile_phylo <- profile_phylo.species_distribution(
    x = as_species_distribution.numeric(x, check_arguments = FALSE),
    tree = tree,
    orders = orders,
    normalize = normalize,
    estimator = estimator,
    level = level,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    gamma = FALSE,
    n_simulations = n_simulations,
    alpha = alpha,
    bootstrap = bootstrap,
    show_progress = show_progress,
    check_arguments = FALSE
  )

  if (as_numeric) {
    return(the_profile_phylo$diversity)
  } else {
    return(the_profile_phylo)
  }
}


#' @rdname profile_phylo
#'
#' @export
profile_phylo.species_distribution <- function(
    x,
    tree,
    orders = seq(from = 0, to = 2, by = 0.1),
    normalize = TRUE,
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger",
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak",
                  "naive"),
    level = NULL,
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05,
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    gamma = FALSE,
    n_simulations = 0,
    alpha = 0.05,
    bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
    show_progress = TRUE,
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
  bootstrap <- match.arg(bootstrap)

  # Calculate abundances along the tree, that are a list of matrices
  the_phylo_abd <- phylo_abd(abundances = x, tree = tree)

  # Prepare arrays to store entropy (3 dimensions: x, y, z)
  # and simulated entropies (4 dimensions : x, y, z, t)
  # x are tree intervals, y are communities, z are orders,
  # t are simulations.
  # Add an array to store simulation envelopes, where
  # t are quantiles of simulations, inf and sup.
  if (gamma) {
    n_communities <- 1
  } else {
    n_communities <- nrow(x)
  }
  ent_phylo_abd <- array(
    dim = c(length(tree$intervals), nrow(x), length(orders))
  )
  if (n_simulations > 0) {
    ent_phylo_sim <- array(
      dim = c(length(tree$intervals), nrow(x), length(orders), n_simulations)
    )
    ent_phylo_envelope <- array(
      dim = c(length(tree$intervals), nrow(x), length(orders), 2)
    )
  }

  # Prepare the progress bar
  if (show_progress & interactive()) {
    cli::cli_progress_bar(
      "Computing phyloentropy",
      total = length(the_phylo_abd) * length(orders)
    )
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
            bootstrap = bootstrap,
            check_arguments = FALSE
          )
        }
      )
      # Prepare an array to store simulated abd. Rows are species, columns are
      # communities (same structure as groups of the_phylo_abd), z are simulations
      # Max number of species in simulated communities
      sp_sim <- max(vapply(comm_sim.list, dim, c(0L,0L))[2,])
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
    for (z_order in seq_along(orders)) {
      # Actual data
      ent_phylo_abd[x_interval, , z_order] <- ent_tsallis.species_distribution(
        x = as_abundances.numeric(t(the_phylo_abd[[x_interval]])),
        q = orders[z_order],
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
      )$entropy

      for (t_simulation in seq_len(n_simulations)) {
        # Entropy of simulated communities
        ent_phylo_sim[x_interval, , z_order, t_simulation] <- ent_tsallis.species_distribution(
          x = as_abundances.numeric(t(comm_sim[, , t_simulation])),
          q = orders[z_order],
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
        )$entropy
      }
      if (n_simulations > 0) {
        for (y_community in seq_len(n_communities)) {
          # Quantiles, recentered
          ent_phylo_envelope[x_interval, y_community, z_order, ] <- stats::quantile(
            ent_phylo_sim[x_interval, y_community, z_order, ],
            probs = c(alpha / 2, 1 - alpha / 2),
            na.rm = TRUE
          ) - mean(ent_phylo_sim[x_interval, y_community, z_order, ]) +
            ent_phylo_abd[x_interval, y_community, z_order]
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
  the_profile_phylo <- div.tibble(
    ent.matrix = ent_community,
    x = x,
    orders = orders
  )
  if (n_simulations > 0) {
    div_inf <- div.tibble(
      ent.matrix = ent_quantiles[, , 1],
      x = x,
      orders = orders
    )
    div_sup <- div.tibble(
      ent.matrix = ent_quantiles[, , 2],
      x = x,
      orders = orders
    )
    the_profile_phylo <- dplyr::bind_cols(
      the_profile_phylo,
      inf = div_inf$diversity,
      sup = div_sup$diversity
    )
  }

  class(the_profile_phylo) <- c("profile", class(the_profile_phylo))
  return(the_profile_phylo)
}


#' Make a long tibble of diversity with a matrix of entropy
#'
#' Utility for [profile_phylo.species_distribution]
#'
#' @param ent.matrix The matrix of entropies.
#' Rows are communities, columns are orders of diversity.
#' @param x The species distribution.
#' @param orders The orders of diversity
#'
#' @returns A tibble. Columns are "site", "order" and "diversity".
#' @noRd
#'
div.tibble <- function(ent.matrix, x, orders) {

  if (!is.matrix(ent.matrix)) {
    # ent.matrix may be a numeric vector (single community / min and max)
    ent.matrix <- t(as.matrix(ent.matrix))
  }
  # Make a tibble with site names and entropies.
  # Columns are orders of diversity
  ent.tibble <- tibble::tibble(
    site = x$site,
    data.frame(ent.matrix)
  )
  colnames(ent.tibble)[-1] <- as.character(orders)
  # Make a long tibble with an "order" column
  ent.tibble <- tidyr::pivot_longer(
    ent.tibble,
    cols = !.data$site,
    names_to = "order",
    names_transform = list(order = as.numeric),
    values_to = "entropy"
  )
  # Calculate diversity
  the_div.tibble <- dplyr::mutate(
    ent.tibble,
    diversity = exp_q(.data$entropy, q = .data$order),
    .keep = "all"
  )
  # Eliminate the "entropy" column
  the_div.tibble <- dplyr::select(
    the_div.tibble,
    -.data$entropy
  )
  # Return
  return(the_div.tibble)
}
