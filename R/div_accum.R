#' Diversity Accumulation of a Community
#' 
#' Details TODO. 
#'
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities]
#' @param ... Unused.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#'
#' @return A tibble with the site names, the estimators used and the accumulated entropy
#' or diversity at each level of sampling effort.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' autoplot(div_accum(paracou_6_abd))
#' 
#' @name div_accum
NULL


#' @rdname div_accum
#'
#' @export
ent_accum <- function(x, ...) {
  UseMethod("ent_accum")
}


#' @rdname div_accum
#'
#' @param q The order of diversity.
#' @param levels The levels, i.e. the sample sizes of interpolation or extrapolation: a vector of integer values.
#' @param probability_estimator One of the estimators of a probability distribution: 
#' "naive" (the default value), or "Chao2013", "Chao2015", "ChaoShen" to estimate
#' the probabilities of the observed species in the asymptotic distribution.
#' @param unveiling One of the possible unveiling methods to estimate the probabilities 
#' of the unobserved species: "none" (default, no species is added), "uniform" 
#' (all unobserved species have the same probability) or "geometric" (the 
#' unobserved species distribution is geometric).
#' @param richness_estimator An estimator of richness to evaluate the total number of species,
#' see [div_richness].
#' @param jack_alpha The risk level, 5% by default, used to optimize the jackknife order.
#' @param jack_max The highest jackknife order allowed. Default is 10. 
#' @param n_simulations The number of simulations used to estimate the confidence envelope of the accumulation.
#' @param alpha The risk level, 5% by default, of the confidence envelope of the accumulation.
#' @param show_progress If TRUE, a progress bar is shown during long computations. 
#' 
#' @export
ent_accum.numeric <- function(
    x, 
    q = 0,
    levels = seq_len(sum(x)), 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  if (!is_integer_values(x)) {
    warning(
      "Integer abundance values are required to estimate community probabilities. Abundances have been rounded."
    )
    x <- round(x)
  }
  
  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)
  # Probabilities
  prob <- abd / sample_size
  # Number of observed species
  s_obs <- length(abd)
  
  # Prepare the vector of results
  ent_level <- numeric(length(levels))
  ent_estimator <- character(length(levels))
  # Prepare the progress bar
  pgb <- utils::txtProgressBar(min = 0, max = length(levels))
  # i must be initialized if the accumulation contains extrapolation only
  i <- 0
  
  # Interpolation ----
  levels_interp <- levels[levels < sample_size]
  # Calculate entropy at each level
  for(level in levels_interp) {
    # Calculate Entropy
    i <- which(levels == level)
    ent_level[i] <- ent_tsallis.numeric(
      abd, 
      q = q, 
      level = level, 
      check_arguments = FALSE
    )
    if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
  }
  ent_estimator[levels_interp] <- "Interpolation"
  
  # level == Sample Size ----
  if (any(levels == sample_size)) {
    i <- which(levels == sample_size)
    ent_level[i] <- ent_tsallis.numeric(
      prob, 
      q = q,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    ent_estimator[i] <- "Sample"
    if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
  }
  
  # Extrapolation ----
  # Don't use Tsallis for speed: probabily extrapolation should be run once only.
  levels_extrap <- levels[levels > sample_size]
  prob_unv <- NULL
  if (length(levels_extrap) > 0) {
    # Unveil the full distribution that rarefies to the observed entropy (or other options)
    prob_unv <- probabilities.numeric(
      abd,
      estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator, 
      jack_alpha  = 0.05, 
      jack_max = 10,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    if (q == 0) {
      ## Richness ----
      s_1 <-  sum(abd == 1)
      if (s_1) {
        # Estimate the number of unobserved species
        s_obs <- sum(abd > 0)
        s_0 <- length(prob_unv) - s_obs
        # Extrapolate richness (the vector is levels_extrap)
        ent_level[(i + 1):length(levels)] <- s_obs + 
          s_0 * (1 - (1 - s_1 / (sample_size * s_0 + s_1))^(levels_extrap - sample_size)) - 1
      } else {
        # No singleton
        ent_level[(i + 1):length(levels)] <- s_obs - 1
      }
      ent_estimator[(i + 1):length(levels)] <- richness_estimator
      if(show_progress & interactive()) utils::setTxtProgressBar(pgb, length(levels))
    } else {
      ## Shannon ----
      if (q == 1) {
        # Estimate the asymptotic entropy
        ent_est <- -sum(prob_unv * log(prob_unv))
        # Estimate observed entropy
        ent_obs <- -sum(prob * log(prob))
        # Interpolation (the vector is levels_extrap)
        ent_level[(i+1):length(levels)] <- sample_size / levels_extrap * ent_obs + 
          (levels_extrap - sample_size) / levels_extrap * ent_est
        ent_estimator[(i + 1):length(levels)] <- richness_estimator
        if(show_progress & interactive()) utils::setTxtProgressBar(pgb, length(levels))
      } else {
        ## Simpson ----
        if (q == 2) {
          # Exit if abd contains no or a single species
          if (s_obs < 2) {
            if (s_obs == 0) {
              ent_level[(i + 1):length(levels)] <- NA
            } else {
              ent_level[(i + 1):length(levels)] <- 0
            }
          } else {
            # Valid extrapolation (the vector is levels_extrap)
            ent_level[(i + 1):length(levels)] <- 1 - 1 / levels_extrap - 
              (1 - 1 / levels_extrap) * sum(abd * (abd - 1)) / sample_size / (sample_size - 1)
          }
          ent_estimator[(i + 1):length(levels)] <- "Chao2014"
          if(show_progress & interactive()) utils::setTxtProgressBar(pgb, length(levels))
        } else {
          # General case: q is not 0, 1 or 2 ----
          for(level in levels_extrap) {
            # Abundance frequence count at level (Chao et al., 2014, eq. 5)
            s_nu <- vapply(
              seq_len(level), 
              function(nu) {
                sum(
                  exp(
                    lchoose(level, nu) + nu * log(prob_unv) + (level - nu) * log(1 - prob_unv)
                  )
                )
              }, 
              FUN.VALUE=0.0
            )
            # Estimate entropy (Chao et al., 2014, eq. 6)
            i <- which(levels == level)
            ent_level[i] <- (sum((seq_len(level) / level)^q * s_nu) - 1) / (1 - q)
            ent_estimator[(i + 1):length(levels)] <- richness_estimator
            if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
          }
        }
      }
    }
  }

  # Simulations ----
  # Generate distributions from the unveiled probabilities
  if (n_simulations > 0 & (probability_estimator == "naive" | unveiling == "none")) {
    warning("Accumulation confidence interval can't be estimated without unveiling the asymptotic distribution. 'probability_estimator' can't be 'naive' and 'unveiling' can't be 'none'.")
    n_simulations <- 0
  }
  if (n_simulations > 0) {
    # Prepare the result matrix
    ent_sim_quantiles <- matrix(0, nrow = length(levels), ncol = 2)
    if (is.null(prob_unv)) {
      # Unveil the full distribution if not done before
      prob_unv <- probabilities.numeric(
        abd,
        estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator, 
        jack_alpha  = 0.05, 
        jack_max = 10,
        as_numeric = TRUE,
        check_arguments = FALSE
      )
    }
    for(level in levels) {
      # Generate simulated communities at each level
      communities <- stats::rmultinom(n_simulations, size = level, prob = prob_unv)
      # Probabilities
      communities <- communities/level
      # Calculate entropy
      ent_sim <- apply(
        communities, 
        MARGIN = 2, 
        FUN = ent_tsallis.numeric, 
        q = q,
        as_numeric = TRUE,      
        check_arguments = FALSE
      )
      i <- which(levels == level)
      # Store quantiles
      ent_sim_quantiles[i, ] <- stats::quantile(
        ent_sim, 
        probs = c(alpha / 2, 1 - alpha / 2)
      )
      if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
    }
  }
  close(pgb)

  # Format the result ----
  if (n_simulations > 0) {
    the_ent_accum <- tibble::tibble(
      level = levels,
      estimator = ent_estimator,
      entropy = ent_level,
      inf = ent_sim_quantiles[, 1],
      sup = ent_sim_quantiles[, 2]
    )
  } else {
    the_ent_accum <- tibble::tibble(
      level = levels,
      estimator = ent_estimator,
      entropy = ent_level
    )
  }
  class(the_ent_accum) <- c("accumulation", class(the_ent_accum))
  
  return(the_ent_accum)
}


#' @rdname div_accum
#'
#' @export
ent_accum.abundances <- function(
    x,
    q = 0,
    levels = NULL, 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Set levels if needed
  if (is.null(levels)) {
    sample_size <- max(
      rowSums(
        x[, !(colnames(x) %in% c("site", "weight"))]
      )
    )
    levels <- seq_len(sample_size)
  }
  # Apply ent_accum.numeric() to each site
  ent_accum_list <- apply(
    # Eliminate site and weight columns
    x[, !(colnames(x) %in% c("site", "weight"))], 
    # Apply to each row
    MARGIN = 1,
    FUN = ent_accum.numeric,
    # Arguments
    q = q,
    levels = levels, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )
  
  # Add site names if needed
  if ("site" %in% colnames(x)) {
    site_names <- x$site
  } else {
    site_names <- paste("site", seq_len(nrow(x)), sep = "_")
  }
  
  # Make a tibble with site, level and entropy
  the_ent_accum <- tibble::tibble(
    site = rep(site_names, each = length(levels)),
    # Coerce the list returned by apply into a dataframe
    do.call(rbind.data.frame, ent_accum_list)
  )
  class(the_ent_accum) <- c("accumulation", class(the_ent_accum))
  
  return(the_ent_accum)
}


#' @rdname div_accum
#'
#' @export
div_accum <- function(x, ...) {
  UseMethod("div_accum")
}


#' @rdname div_accum
#'
#' @export
div_accum.numeric <- function(
    x, 
    q = 0,
    levels = seq_len(sum(x)), 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Accumulate entropy
  the_div_accum <- ent_accum.numeric(
    x,
    q = q,
    levels = levels, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )
  
  # Calculate diversity
  the_div_accum <- dplyr::mutate(
    the_div_accum,
    diversity = exp_q(.data$entropy, q = q)
  )

  return(the_div_accum)
}


#' @rdname div_accum
#'
#' @export
div_accum.abundances <- function(
    x,
    q = 0,
    levels = NULL, 
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Set levels if needed
  if (is.null(levels)) {
    sample_size <- max(
      rowSums(
        x[, !(colnames(x) %in% c("site", "weight"))]
      )
    )
    levels <- seq_len(sample_size)
  }
  # Apply ent_accum.numeric() to each site
  div_accum_list <- apply(
    # Eliminate site and weight columns
    x[, !(colnames(x) %in% c("site", "weight"))], 
    # Apply to each row
    MARGIN = 1,
    FUN = div_accum.numeric,
    # Arguments
    q = q,
    levels = levels, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )
  
  # Add site names if needed
  if ("site" %in% colnames(x)) {
    site_names <- x$site
  } else {
    site_names <- paste("site", seq_len(nrow(x)), sep = "_")
  }
  
  # Make a tibble with site, level and diversity
  the_div_accum <- tibble::tibble(
    site = rep(site_names, each = length(levels)),
    # Coerce the list returned by apply into a dataframe
    do.call(rbind.data.frame, div_accum_list)
  )
  class(the_div_accum) <- c("accumulation", class(the_div_accum))
  
  return(the_div_accum)
}


#' @rdname div_accum
#'
#' @param object An object of class "accumulation".
#' @param main The main title of the plot.
#' @param xlab The label of the x-axis.
#' @param ylab The label of the y-axis.
#' @param shade_color The color of the shaded confidence envelopes.
#' @param alpha The opacity of the confidence envelopes, between 0 (transparent) and 1 (opaque).
#' @param border_color The color of the border of the confidence envelopes.
#' @param lty The line type of the curves.
#' @param lwd The line width of the curves.
#'
#' @export
autoplot.accumulation <-  function(
    object, 
    ..., 
    main = NULL,
    xlab = "Sample Size",
    ylab = NULL,
    shade_color = "grey75",
    alpha = 0.3,
    border_color = "red",
    lty = ggplot2::GeomLine$default_aes$linetype,
    lwd = ggplot2::GeomLine$default_aes$linewidth){
  
  # Add a site column if needed
  if (!("site" %in% colnames(object))) {
    object <- dplyr::mutate(object, site = "Unique site")
  }

  if ("diversity" %in% colnames(object)) {
    the_plot <- ggplot2::ggplot(
      object, 
      ggplot2::aes(
        x = .data$level, 
        y = .data$diversity,
        col = .data$site
      )
    )
  } else {
    the_plot <- ggplot2::ggplot(
      object, 
      ggplot2::aes(
        x = .data$level, 
        y = .data$entropy,
        color = .data$site
      )
    )
  }
  if ("sup" %in% colnames(object) & "inf" %in% colnames(object)) {
    the_plot <- the_plot +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$inf, 
          ymax = .data$sup
        ), 
        fill = shade_color, 
        alpha = alpha
      ) +
      # Add red lines on borders of polygon
      ggplot2::geom_line(
        ggplot2::aes(y = .data$inf), 
        color = border_color, 
        linetype = 2
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$sup), 
        color = border_color, 
        linetype = 2
      )
  }
  if (is.null(ylab)) {
    if ("diversity" %in% colnames(object)) {
      ylab <- "Diversity"
    } else {
      ylab <- "Entropy"
    }
  }
  the_plot <- the_plot +
    ggplot2::geom_line(linetype = lty, linewidth = lwd) +
    ggplot2::labs(title = main, x = xlab, y = ylab)
  
  # Actual value
  actual <- dplyr::filter(object, .data$estimator == "Sample")
  if (nrow(actual) > 0) {
    if ("diversity" %in% colnames(object)) {
      the_plot <- the_plot +
        ggplot2::geom_hline(
          data = actual, 
          mapping = ggplot2::aes(yintercept = .data$diversity, color = .data$site), 
          linetype = 2
        )
    } else {
      the_plot <- the_plot +
        ggplot2::geom_hline(
          data = actual, 
          mapping = ggplot2::aes(yintercept = .data$entropy, color = .data$site), 
          linetype = 2
        )
    }
    the_plot <- the_plot +
      ggplot2::geom_vline(
        data = actual, 
        mapping = ggplot2::aes(xintercept = .data$level, color = .data$site), 
        linetype = 2
      ) 
  }
  
  return(the_plot)
}
