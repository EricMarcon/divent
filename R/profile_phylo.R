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
#' @return A tibble with the site names, the estimators used and the estimated diversity at each order.
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
#' @param bootstrap The method used to obtain the probabilities to generate 
#' bootstrapped communities from observed abundances. 
#' If "Marcon2012", the probabilities are simply the abundances divided by the total
#' number of individuals \insertCite{Marcon2012a}{divent}. 
#' If "Chao2013" or "Chao2015" (by default), a more sophisticated approach is used 
#' (see [as_probabilities]) following \insertCite{Chao2013;textual}{divent} or 
#' \insertCite{Chao2015;textual}{divent}.
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
    x = as_species_distribution.numeric(x), 
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
  
  # Prepare an array to store diversity (3 dimensions: x, y, z)
  # and simulated entropies (4 dimensions : x, y, z, t)
  # x are tree intervals, y are communities, z are orders, 
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
      dim = c(length(tree$intervals), nrow(x), length(orders), 2)
    )
  }
  
  # Prepare the progress bar
  pgb <- utils::txtProgressBar(
    min = 0, 
    max = length(the_phylo_abd) * n_communities * length(orders)
  )
  
  # Calculate entropy along the tree
  for (x_interval in seq_along(the_phylo_abd)) {
    # Intervals are items of the list
    if (gamma) {
      abd <- metacommunity(
        the_phylo_abd[[x_interval]] # TODO : all arguments
      )
      # TODO : metacommunity.numeric. ! abd must be a column matrix
    } else {
      abd <- the_phylo_abd[[x_interval]]
    }
    for (y_community in seq_len(n_communities)) {
      # Abundances are in columns: abd[, y_community]
      if (n_simulations > 0) {
        # Simulate communities
        comm_sim <- rcommunity(
          n_simulations,
          abd = abd[, y_community],
          bootstrap = bootstrap,
          check_arguments = FALSE
        )
      }
      for (z_order in seq_along(orders)) {
        # Actual data
        ent_phylo_abd[x_interval, y_community, z_order] <- ent_tsallis.numeric(
          x = abd[, y_community],
          q = orders[z_order],
          estimator = estimator,
          level = level, 
          probability_estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator,
          jack_alpha  = jack_alpha, 
          jack_max = jack_max, 
          coverage_estimator = coverage_estimator,
          as_numeric = TRUE,
          check_arguments = FALSE
        )
        if (n_simulations > 0) {
          # Simulated communities
          ent_sim <- ent_tsallis.species_distribution(
            x = comm_sim,
            q = orders[z_order],
            estimator = estimator,
            level = level, 
            probability_estimator = probability_estimator,
            unveiling = unveiling,
            richness_estimator = richness_estimator,
            jack_alpha  = jack_alpha, 
            jack_max = jack_max, 
            coverage_estimator = coverage_estimator,
            check_arguments = FALSE
          )
          # Quantiles, recentered
          ent_phylo_sim[x_interval, y_community, z_order, ] <- stats::quantile(
            ent_sim$entropy, 
            probs = c(alpha / 2, 1 - alpha / 2)
          ) - mean(ent_sim$entropy) + 
            ent_phylo_abd[x_interval, y_community, z_order]
          # Progress bar
        }
      }
      if (show_progress & interactive()) {
        utils::setTxtProgressBar(pgb, utils::getTxtProgressBar(pgb) + 1)
      }
    }
  }
  close(pgb)

  # Average entropy
  # Actual data
  ent_community <- apply(
    ent_phylo_abd,
    MARGIN = 2:3,
    FUN = stats::weighted.mean,
    # Arguments
    w = tree$intervals
  )
  if (!is.matrix(ent_community)) {
    # With a single community, the_entropy is a numeric vector 
    ent_community <- t(as.matrix(ent_community))
  }
  # Simulations
  if (n_simulations > 0) {
    ent_quantiles <- apply(
      ent_phylo_sim,
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


#' Make a long tibble with a matrix
#' 
#' Utility for [profile_phylo.species_distribution]
#'
#' @param ent.matrix The matrix of entropies. 
#' Rows are communities, columns are orders of diversity.
#' @param x The species distribution.
#' @param orders The orders of diversity
#'
#' @return A tibble. Columns are "site", "order" and "diversity".
#' @noRd
#'
div.tibble <- function(ent.matrix, x, orders) {
  
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
    cols = !site,
    names_to = "order",
    names_transform = list(order = as.numeric),
    values_to = "entropy"
  )
  # Calculate diversity
  the_div.tibble <- dplyr::mutate(
    ent.tibble, 
    diversity = exp_q(.data$entropy, q = order),
    .keep = "all"
  )
  # Eliminate the "entropy" column
  the_div.tibble <- dplyr::select(
    the_div.tibble,
    -entropy
  )
  # Return
  return(the_div.tibble)
}
