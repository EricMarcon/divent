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
  
  # Numeric vector, no simulation ----
  if (as_numeric) {
    if (n_simulations > 0) stop ("No simulations are allowed if a numeric vector is expected ('as_numeric = TRUE').")
    the_profile_phylo <- vapply(
      orders,
      FUN = function(q) {
        div_phylo.numeric(
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
          as_numeric = TRUE,
          check_arguments = FALSE
        )
      },
      FUN.VALUE = 0
    )
    return(the_profile_phylo)
  } 
  
  # Regular output, simulations are allowed ----
  the_profile_phylo <- lapply(
    orders,
    FUN = function(q) {
      div_phylo.numeric(
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
    }
  )
  # Make a tibble with the list
  the_profile_phylo <- do.call(rbind.data.frame, the_profile_phylo)
  
  if (n_simulations > 0) {
    # Simulations ----
    if (!is_integer_values(x)) {
      warning("Evaluation of the confidence interval of community profiles requires integer abundances. They have been rounded.")
    }
    abd_int <- round(x)
    # Simulate communities
    communities <- rcommunity(
      n_simulations,
      abd = abd_int,
      bootstrap = bootstrap,
      check_arguments = FALSE
    )
    # Prepare the progress bar
    pgb <- utils::txtProgressBar(min = 0, max = n_simulations)
    # Prepare the result matrix
    profile_phylos <- matrix(0, nrow = n_simulations, ncol = length(orders))
    # Loops are required for the progress bar
    for (i in seq_len(n_simulations)) {
      # Parallelize. Do not allow more forks.
      profiles_list <- parallel::mclapply(
        orders, 
        FUN = function(q) {
          div_phylo.numeric(
            communities[i, !colnames(communities) %in% non_species_columns], 
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
            as_numeric = TRUE,
            check_arguments = FALSE
          )
        },
        mc.allow.recursive = FALSE
      )
      profile_phylos[i, ] <- simplify2array(profiles_list)
      if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
    }
    close(pgb)
    # Recenter simulated values
    div_means <- apply(profile_phylos, 2, mean)
    profile_phylos <- t(
      t(profile_phylos) - div_means + the_profile_phylo$diversity
    )
    # Quantiles
    div_quantiles <- apply(
      profile_phylos, 
      MARGIN = 2, 
      FUN = stats::quantile,
      probs = c(alpha / 2, 1 - alpha / 2)
    )
    # Format the result 
    the_profile_phylo <- tibble::tibble(
      the_profile_phylo,
      inf = div_quantiles[1, ],
      sup = div_quantiles[2, ]
    )
  }
  class(the_profile_phylo) <- c("profile", class(the_profile_phylo))
  
  return(the_profile_phylo)
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
  
  if (gamma) {
    the_profile_phylo <- profile_phylo.numeric(
      metacommunity(x, as_numeric = TRUE, check_arguments = FALSE),
      # Arguments
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
      n_simulations = n_simulations,
      alpha = alpha,
      bootstrap = bootstrap,
      show_progress = show_progress,
      check_arguments = FALSE
    )
  } else {
    # Apply profile_phylo.numeric() to each site
    profile_phylo_list <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns], 
      # Apply to each row
      MARGIN = 1,
      FUN = profile_phylo.numeric,
      # Arguments
      tree = tree,
      orders = orders,
      normalize = normalize,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max, 
      coverage_estimator = coverage_estimator,
      as_numeric = FALSE,
      n_simulations = n_simulations,
      alpha = alpha,
      bootstrap = bootstrap,
      show_progress = show_progress,
      check_arguments = FALSE
    )
    # Make a tibble with sites and profiles
    the_profile_phylo <- tibble::tibble(
      site = rep(x$site, each = length(orders)),
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, profile_phylo_list)
    )
  }
  class(the_profile_phylo) <- c("profile", class(the_profile_phylo))
  
  return(the_profile_phylo)
}
