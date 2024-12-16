#' Similarity-Based Diversity Profile of a Community
#'
#' Calculate the diversity profile of a community, i.e. its similarity-based diversity
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
#' # Similarity matrix
#' Z <- fun_similarity(paracou_6_fundist)
#' # Profile
#' profile_similarity(paracou_6_abd, similarities = Z, q = 2)
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity at each order.
#' This is an object of class "profile" that can be plotted.
#'
#' @references
#' \insertAllCited{}
#'
#' @name profile_similarity
NULL


#' @rdname profile_similarity
#'
#' @export
profile_similarity <- function(
    x,
    similarities,
    orders = seq(from = 0, to = 2, by = 0.1),
    ...) {
  UseMethod("profile_similarity")
}


#' @rdname profile_similarity
#'
#' @param orders The orders of diversity used to build the profile.
#' @param estimator An estimator of entropy.
#' @param n_simulations The number of simulations used to estimate the confidence envelope of the profile.
#' @param alpha The risk level, 5% by default, of the confidence envelope of the profile.
#'
#' @export
profile_similarity.numeric <- function(
    x,
    similarities = diag(length(x)),
    orders = seq(from = 0, to = 2, by = 0.1),
    estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang",
                  "UnveilC", "UnveiliC", "naive"),
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

  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  coverage_estimator <- match.arg(coverage_estimator)
  bootstrap <- match.arg(bootstrap)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    similarities <- checked_matrix(similarities, x)
  }

  # Numeric vector, no simulation ----
  if (as_numeric) {
    if (n_simulations > 0) {
      stop("No simulations are allowed if a numeric vector is expected ('as_numeric = TRUE').")
    }
    the_profile_similarity <- vapply(
      orders,
      FUN = function(q) {
        div_similarity.numeric(
          x,
          similarities = similarities,
          q = q,
          estimator = estimator,
          probability_estimator = probability_estimator,
          unveiling = unveiling,
          jack_alpha  = jack_alpha,
          jack_max = jack_max,
          coverage_estimator = coverage_estimator,
          sample_coverage = sample_coverage,
          as_numeric = TRUE,
          check_arguments = FALSE
        )
      },
      FUN.VALUE = 0
    )
    return(the_profile_similarity)
  }

  # Regular output, simulations are allowed ----
  the_profile_similarity <- lapply(
    orders,
    FUN = function(q) {
      div_similarity.numeric(
        x,
        similarities = similarities,
        q = q,
        estimator = estimator,
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        jack_alpha  = jack_alpha,
        jack_max = jack_max,
        coverage_estimator = coverage_estimator,
        as_numeric = FALSE,
        check_arguments = FALSE
      )
    }
  )
  # Make a tibble with the list
  the_profile_similarity <- do.call(rbind.data.frame, the_profile_similarity)

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
    if (show_progress & interactive()) {
      cli::cli_progress_bar("Running simulations", total = n_simulations)
    }
    # Prepare the result matrix
    profile_similarities <- matrix(0, nrow = n_simulations, ncol = length(orders))
    # Loops are required for the progress bar
    for (i in seq_len(n_simulations)) {
      # Parallelize. Do not allow more forks.
      profiles_list <- parallel::mclapply(
        orders,
        FUN = function(q) {
          div_similarity.numeric(
            communities[i, !colnames(communities) %in% non_species_columns],
            similarities = similarities,
            q = q,
            estimator = estimator,
            probability_estimator = probability_estimator,
            unveiling = unveiling,
            jack_alpha  = jack_alpha,
            jack_max = jack_max,
            coverage_estimator = coverage_estimator,
            sample_coverage = sample_coverage,
            as_numeric = TRUE,
            check_arguments = FALSE
          )
        },
        mc.allow.recursive = FALSE
      )
      profile_similarities[i, ] <- simplify2array(profiles_list)
      if (show_progress & interactive()) cli::cli_progress_update()
    }
    # Recenter simulated values
    div_means <- apply(profile_similarities, 2, mean)
    profile_similarities <- t(
      t(profile_similarities) - div_means + the_profile_similarity$diversity
    )
    # Quantiles
    div_quantiles <- apply(
      profile_similarities,
      MARGIN = 2,
      FUN = stats::quantile,
      probs = c(alpha / 2, 1 - alpha / 2)
    )
    # Format the result
    the_profile_similarity <- tibble::tibble(
      the_profile_similarity,
      inf = div_quantiles[1, ],
      sup = div_quantiles[2, ]
    )
  }
  class(the_profile_similarity) <- c("profile", class(the_profile_similarity))

  return(the_profile_similarity)
}


#' @rdname profile_similarity
#'
#' @export
profile_similarity.species_distribution <- function(
    x,
    similarities = diag(sum(!colnames(x) %in% non_species_columns)),
    orders = seq(from = 0, to = 2, by = 0.1),
    estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang",
                  "UnveilC", "UnveiliC", "naive"),
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
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

  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  coverage_estimator <- match.arg(coverage_estimator)
  bootstrap <- match.arg(bootstrap)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
    similarities <- checked_matrix(similarities, x)
  }

  if (gamma) {
    the_profile_similarity <- profile_similarity.numeric(
      metacommunity.abundances(
        x = x,
        as_numeric = TRUE,
        check_arguments = FALSE
      ),
      # Arguments
      similarities = similarities,
      orders = orders,
      estimator = estimator,
      probability_estimator = probability_estimator,
      unveiling = unveiling,
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
    # Apply profile_similarity.numeric() to each site
    profile_similarity_list <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns],
      # Apply to each row
      MARGIN = 1,
      FUN = profile_similarity.numeric,
      # Arguments
      similarities = similarities,
      orders = orders,
      probability_estimator = probability_estimator,
      unveiling = unveiling,
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
    the_profile_similarity <- tibble::tibble(
      site = rep(x$site, each = length(orders)),
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, profile_similarity_list)
    )
  }
  class(the_profile_similarity) <- c("profile", class(the_profile_similarity))

  return(the_profile_similarity)
}
