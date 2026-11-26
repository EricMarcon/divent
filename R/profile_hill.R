#' Diversity Profile of a Community
#'
#' Calculate the diversity profile of a community, i.e. diversity (Hill numbers)
#' against its order.
#'
#' A bootstrap confidence interval can be produced by simulating communities
#' (their number is `n_simulations`) with [rcommunity] and calculating their profiles.
#' Simulating communities implies a downward bias in the estimation:
#' rare species of the actual community may have abundance zero in simulated communities.
#' Simulated diversity values are recentered so that their mean is that of the actual community.
#'
#' Simulations can be run in parallel thanks to package [doFuture].
#' A parallel plan must be chosen with function [future::plan()].
#' To use all cores but one on a personal computer, run:
#'
#' ```
#' library("doFuture")
#' plan(multisession, workers = availableCores(omit = 1))
#' ```
#' If argument `show_progress` is `TRUE`, a progress bar can be displayed
#' by package [progressr].
#' Since all progress bars rely on package [cli] in *divent*, the "cli" handler
#' is a good choice.
#'
#' ```
#' library("progressr")
#' handlers(global = TRUE)
#' handlers("cli")
#' ```
#' Those code chunks must be run before `profile_hill()`
#' with large datasets.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @examples
#' autoplot(profile_hill(paracou_6_abd))
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity at each order.
#' This is an object of class "profile" that can be plotted.
#'
#' If `as_numeric` is `TRUE`, then the value is a matrix whose first column is the order
#' and other columns are the diversity of sites.
#'
#' @references
#' \insertAllCited{}
#'
#' @name profile_hill
NULL


#' @rdname profile_hill
#'
#' @export
profile_hill <- function(
  x,
  orders = seq(from = 0, to = 2, by = 0.1),
  ...
) {
  UseMethod("profile_hill")
}


#' @rdname profile_hill
#'
#' @param orders The orders of diversity used to build the profile.
#' @param estimator An estimator of entropy.
#' @param n_simulations The number of simulations used to estimate the confidence envelope of the profile.
#' @param alpha The risk level, 5% by default, of the confidence envelope of the profile.
#'
#' @export
profile_hill.numeric <- function(
  x,
  orders = seq(from = 0, to = 2, by = 0.1),
  estimator = c(
    "UnveilJ",
    "ChaoJost",
    "ChaoShen",
    "GenCov",
    "Grassberger",
    "Holste",
    "Marcon",
    "UnveilC",
    "UnveiliC",
    "ZhangGrabchak",
    "naive"
  ),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q_threshold = 10,
  sample_coverage = NULL,
  as_numeric = FALSE,
  n_simulations = 0,
  alpha = 0.05,
  bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
) {
  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
  bootstrap <- match.arg(bootstrap)
  if (check_arguments) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }
  if (as_numeric && n_simulations > 0) {
    cli::cli_abort(
      c(
        "No simulations are allowed if a numeric vector is expected",
        "i" = "Change argument for {.code as_numeric = FALSE}"
      )
    )
  }

  # Compute the profile ----
  the_profile_hill <- lapply(
    orders,
    FUN = function(q) {
      div_hill(
        x,
        q = q,
        estimator = estimator,
        level = level,
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha = jack_alpha,
        jack_max = jack_max,
        coverage_estimator = coverage_estimator,
        q_threshold = q_threshold,
        as_numeric = FALSE,
        check_arguments = FALSE
      )
    }
  )
  # Make a tibble with the list
  the_profile_hill <- dplyr::bind_rows(the_profile_hill)

  if (n_simulations > 0) {
    # Simulations ----
    if (!is_integer_values(x)) {
      cli::cli_alert_warning(
        paste(
          "Evaluation of the confidence interval of community profiles requires",
          "integer abundances."
        )
      )
      cli::cli_alert("They have been rounded.")
    }
    abd_int <- round(x)
    # Simulate communities
    communities <- rcommunity(
      n_simulations,
      abd = abd_int,
      bootstrap = bootstrap,
      check_arguments = FALSE
    )

    # Define the progressor
    if (show_progress && interactive()) {
      progress <- progressr::progressor(steps = n_simulations)
    }

    # Run simulations in parallel
    profile_hills <- foreach::foreach(
      i = seq_len(n_simulations),
      .combine = rbind
    ) %dofuture%
      {
        # Calculate diversity at all orders
        profiles_list <- lapply(
          orders,
          FUN = function(q) {
            div_hill(
              communities[i, ],
              q = q,
              estimator = estimator,
              level = level,
              probability_estimator = probability_estimator,
              unveiling = unveiling,
              richness_estimator = richness_estimator,
              jack_alpha = jack_alpha,
              jack_max = jack_max,
              coverage_estimator = coverage_estimator,
              q_threshold = q_threshold,
              as_numeric = TRUE,
              check_arguments = FALSE
            )
          }
        )
        # Progressor
        if (show_progress && interactive()) {
          progress()
        }
        # Return a vector
        simplify2array(profiles_list)
      }
    # Recenter simulated values
    div_means <- colMeans(profile_hills)
    profile_hills <- t(
      t(profile_hills) - div_means + the_profile_hill$diversity
    )
    # Quantiles
    div_quantiles <- apply(
      profile_hills,
      MARGIN = 2,
      FUN = stats::quantile,
      probs = c(alpha / 2, 1 - alpha / 2)
    )
    # Format the result
    the_profile_hill <- tibble::tibble(
      the_profile_hill,
      inf = div_quantiles[1, ],
      sup = div_quantiles[2, ]
    )
  }

  if (as_numeric) {
    return(the_profile_hill$diversity)
  } else {
    class(the_profile_hill) <- c("profile", class(the_profile_hill))
    return(the_profile_hill)
  }
}


#' @rdname profile_hill
#'
#' @export
profile_hill.species_distribution <- function(
  x,
  orders = seq(from = 0, to = 2, by = 0.1),
  estimator = c(
    "UnveilJ",
    "ChaoJost",
    "ChaoShen",
    "GenCov",
    "Grassberger",
    "Holste",
    "Marcon",
    "UnveilC",
    "UnveiliC",
    "ZhangGrabchak",
    "naive"
  ),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q_threshold = 10,
  gamma = FALSE,
  as_numeric = FALSE,
  n_simulations = 0,
  alpha = 0.05,
  bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
) {
  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
  bootstrap <- match.arg(bootstrap)
  if (check_arguments) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }
  if (as_numeric && n_simulations > 0) {
    cli::cli_abort(
      c(
        "No simulations are allowed if a numeric vector is expected",
        "i" = "Change argument for {.code as_numeric = FALSE}"
      )
    )
  }

  if (gamma) {
    the_profile_hill <- profile_hill(
      metacommunity(
        x = x,
        as_numeric = TRUE,
        check_arguments = FALSE
      ),
      # Arguments
      q = q,
      estimator = estimator,
      level = level,
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha = jack_alpha,
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      q_threshold = q_threshold,
      as_numeric = FALSE,
      n_simulations = n_simulations,
      alpha = alpha,
      bootstrap = bootstrap,
      show_progress = show_progress,
      check_arguments = FALSE
    )
  } else {
    # Apply profile_hill.numeric() to each site
    profile_hill_list <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns],
      # Apply to each row
      MARGIN = 1,
      FUN = profile_hill,
      # Arguments
      orders = orders,
      estimator = estimator,
      level = level,
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha = jack_alpha,
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      q_threshold = q_threshold,
      as_numeric = FALSE,
      n_simulations = n_simulations,
      alpha = alpha,
      bootstrap = bootstrap,
      show_progress = show_progress,
      check_arguments = FALSE
    )
    # Make a tibble with sites and profiles
    the_profile_hill <- tibble::tibble(
      site = rep(x$site, each = length(orders)),
      # Coerce the list returned
      dplyr::bind_rows(profile_hill_list)
    )
  }
  if (as_numeric) {
    # Make a matrix: first col is order, other cols are diversity of sites
    return(as.matrix(pivot_wider(
      the_profile_hill,
      names_from = site,
      values_from = diversity
    )[, -1]))
  } else {
    class(the_profile_hill) <- c("profile", class(the_profile_hill))
    return(the_profile_hill)
  }
}
