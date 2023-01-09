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
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @examples
#' profile_hill(paracou_6_abd)
#'
#' @return A tibble with the site names, the estimators used and the estimated diversity at each order.
#' This is an object of class "profile" that can be plotted.
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
    ...) {
  UseMethod("profile_hill")
}


#' @rdname profile_hill
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
profile_hill.numeric <- function(
    x, 
    orders = seq(from = 0, to = 2, by = 0.1), 
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
    the_profile_hill <- vapply(
      orders,
      FUN = function(q) {
        div_hill.numeric(
          x,
          q = q,
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
    return(the_profile_hill)
  } 
  
  # Regular output, simulations are allowed ----
  the_profile_hill <- lapply(
    orders,
    FUN = function(q) {
      div_hill.numeric(
        x,
        q = q,
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
  the_profile_hill <- do.call(rbind.data.frame, the_profile_hill)
  
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
    profile_hills <- matrix(0, nrow = n_simulations, ncol = length(orders))
    # Loops are required for the progress bar
    for (i in seq_len(n_simulations)) {
      # Parallelize. Do not allow more forks.
      profiles_list <- parallel::mclapply(
        orders, 
        FUN = function(q) {
          div_hill.numeric(
            communities[i, !colnames(communities) %in% non_species_columns], 
            q = q,
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
      profile_hills[i, ] <- simplify2array(profiles_list)
      if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
    }
    close(pgb)
    # Recenter simulated values
    div_means <- apply(profile_hills, 2, mean)
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
  class(the_profile_hill) <- c("profile", class(the_profile_hill))
  
  return(the_profile_hill)
}


#' @rdname profile_hill
#'
#' @export
profile_hill.species_distribution <- function(
    x, 
    orders = seq(from = 0, to = 2, by = 0.1), 
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
  }
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)
  bootstrap <- match.arg(bootstrap)

  if (gamma) {
    the_profile_hill <- profile_hill.numeric(
      metacommunity.abundances(
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
    # Apply profile_hill.numeric() to each site
    profile_hill_list <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns], 
      # Apply to each row
      MARGIN = 1,
      FUN = profile_hill.numeric,
      # Arguments
      orders = orders,
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
    # Make a tibble with sites and profiles
    the_profile_hill <- tibble::tibble(
      site = rep(x$site, each = length(orders)),
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, profile_hill_list)
    )
  }
  class(the_profile_hill) <- c("profile", class(the_profile_hill))
  
  return(the_profile_hill)
}


#' @rdname profile_hill
#'
#' @param object An object of class "profile".
#' @param main The main title of the plot.
#' @param xlab The label of the x-axis.
#' @param ylab The label of the y-axis.
#' @param shade_color The color of the shaded confidence envelopes.
#' @param alpha The opacity of the confidence envelopes, between 0 (transparent) and 1 (opaque).
#' @param lty The line type of the curves.
#' @param lwd The line width of the curves.
#'
#' @export
autoplot.profile <-  function(
    object, 
    ..., 
    main = NULL,
    xlab = "Order of Diversity",
    ylab = "Diversity",
    shade_color = "grey75",
    alpha = 0.3,
    lty = ggplot2::GeomLine$default_aes$linetype,
    lwd = ggplot2::GeomLine$default_aes$linewidth){
  
  # Add a site column if needed
  if (!"site" %in% colnames(object)) {
    object <- dplyr::mutate(object, site = "Unique site")
  }
  
  # Build the plot
  the_plot <- ggplot2::ggplot(
    object, 
    ggplot2::aes(
      x = .data$order, 
      y = .data$diversity,
      col = .data$site
    )
  )

  # Confidence envelope
  if ("sup" %in% colnames(object) & "inf" %in% colnames(object)) {
  the_plot <- the_plot +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$inf, 
        ymax = .data$sup,
        lty = .data$site
      ), 
      fill = shade_color, 
      alpha = alpha,
      col = shade_color
    )
  }

  # Profiles
  the_plot <- the_plot +
    ggplot2::geom_line(linetype = lty, linewidth = lwd) +
    ggplot2::labs(title = main, x = xlab, y = ylab)
  
  return(the_plot)
}
