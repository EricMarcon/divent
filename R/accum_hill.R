#' Diversity Accumulation of a Community
#'
#' Diversity and Entropy Accumulation Curves represent the accumulation of
#' entropy and diversity with respect to the sample size.
#'
#' `accum_hill()` or `accum_tsallis()` estimate the diversity or entropy accumulation
#' curve of a distribution.
#' See [ent_tsallis] for details about the computation of entropy at each level
#' of interpolation and extrapolation.
#'
#' In accumulation curves, extrapolation is done by estimating the asymptotic
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
#' # Paracou 6 subplot 1
#' autoplot(accum_hill(paracou_6_abd[1, ]))
#'
#' @name accum_hill
NULL


#' @rdname accum_hill
#'
#' @export
accum_tsallis <- function(x, ...) {
  UseMethod("accum_tsallis")
}


#' @rdname accum_hill
#'
#' @param levels The levels, i.e. the sample sizes of interpolation or
#' extrapolation: a vector of integer values.
#'
#' @export
accum_tsallis.numeric <- function(
    x,
    q = 0,
    levels = seq_len(sum(x)),
    probability_estimator = c("Chao2015", "Chao2013","ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha = 0.05,
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    n_simulations = 0,
    alpha = 0.05,
    show_progress = TRUE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)

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
  if (show_progress & interactive()) {
    cli::cli_progress_bar("Running estimations", total = length(levels))
  }
  # i must be initialized if the accumulation contains extrapolation only
  i <- 0

  # Interpolation ----
  levels_interp <- levels[levels < sample_size]
  # Calculate entropy at each level
  for (level in levels_interp) {
    # Calculate Entropy
    i <- which(levels == level)
    ent_level[i] <- ent_tsallis.numeric(
      abd,
      q = q,
      level = level,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    ent_estimator[i] <- "Interpolation"
    if (show_progress & interactive()) cli::cli_progress_update(set = i)
  }

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
    if (show_progress & interactive()) cli::cli_progress_update(set = i)
  }

  # Extrapolation ----
  # Don't use ent_tsallis for speed: probability extrapolation should be run once only.
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
      coverage_estimator = coverage_estimator,
      q = q,
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
      if (show_progress & interactive()) {
        cli::cli_progress_bar("Running estimations", total = length(levels))
      }
    } else {
      ## Shannon ----
      if (q == 1) {
        # Estimate the asymptotic entropy
        ent_est <- -sum(prob_unv * log(prob_unv))
        # Estimate observed entropy
        ent_obs <- -sum(prob * log(prob))
        # Interpolation (the vector is levels_extrap)
        ent_level[(i + 1):length(levels)] <- sample_size / levels_extrap * ent_obs +
          (levels_extrap - sample_size) / levels_extrap * ent_est
        ent_estimator[(i + 1):length(levels)] <- richness_estimator
        if (show_progress & interactive()) {
          cli::cli_progress_bar("Running estimations", total = length(levels))
        }
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
          if (show_progress & interactive()) {
            cli::cli_progress_bar("Running estimations", total = length(levels))
          }
        } else {
          # General case: q is not 0, 1 or 2 ----
          for (level in levels_extrap) {
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
              FUN.VALUE = 0.0
            )
            # Estimate entropy (Chao et al., 2014, eq. 6)
            i <- which(levels == level)
            ent_level[i] <- (sum((seq_len(level) / level)^q * s_nu) - 1) / (1 - q)
            ent_estimator[i] <- richness_estimator
            if (show_progress & interactive()) cli::cli_progress_update(set = i)
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
    if (show_progress & interactive()) {
      cli::cli_progress_bar("Running simulations", total = length(levels))
    }
    if (is.null(prob_unv)) {
      # Unveil the full distribution if not done before
      prob_unv <- probabilities.numeric(
        abd,
        estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = 0.05,
        jack_max = 10,
        coverage_estimator = coverage_estimator,
        q = q,
        as_numeric = TRUE,
        check_arguments = FALSE
      )
    }
    for (level in levels) {
      # Generate simulated communities at each level
      communities <- stats::rmultinom(n_simulations, size = level, prob = prob_unv)
      # Probabilities
      communities <- communities / level
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
      if (show_progress & interactive()) cli::cli_progress_update(set = i)
    }
  }

  # Format the result ----
  if (n_simulations > 0) {
    the_accum_tsallis <- tibble::tibble(
      level = levels,
      estimator = ent_estimator,
      entropy = ent_level,
      inf = ent_sim_quantiles[, 1],
      sup = ent_sim_quantiles[, 2]
    )
  } else {
    the_accum_tsallis <- tibble::tibble(
      level = levels,
      estimator = ent_estimator,
      entropy = ent_level
    )
  }
  class(the_accum_tsallis) <- c("accumulation", class(the_accum_tsallis))

  return(the_accum_tsallis)
}


#' @rdname accum_hill
#'
#' @export
accum_tsallis.abundances <- function(
    x,
    q = 0,
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

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)

  # Set levels if needed
  if (is.null(levels)) {
    sample_size <- max(
      rowSums(
        x[, !colnames(x) %in% non_species_columns]
      )
    )
    levels <- seq_len(sample_size)
  }
  # Apply accum_tsallis.numeric() to each site
  accum_tsallis_list <- apply(
    # Eliminate site and weight columns
    x[, !colnames(x) %in% non_species_columns],
    # Apply to each row
    MARGIN = 1,
    FUN = accum_tsallis.numeric,
    # Arguments
    q = q,
    levels = levels,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
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
  the_accum_tsallis <- tibble::tibble(
    site = rep(site_names, each = length(levels)),
    # Coerce the list returned by apply into a dataframe
    do.call(rbind.data.frame, accum_tsallis_list)
  )
  class(the_accum_tsallis) <- c("accumulation", class(the_accum_tsallis))

  return(the_accum_tsallis)
}


#' @rdname accum_hill
#'
#' @export
accum_hill <- function(x, ...) {
  UseMethod("accum_hill")
}


#' @rdname accum_hill
#'
#' @export
accum_hill.numeric <- function(
    x,
    q = 0,
    levels = seq_len(sum(x)),
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

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)

  # Accumulate entropy
  the_accum_tsallis <- accum_tsallis.numeric(
    x,
    q = q,
    levels = levels,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    n_simulations = n_simulations,
    alpha = alpha,
    show_progress = show_progress,
    check_arguments = FALSE
  )

  # Calculate diversity
  the_accum_hill <- dplyr::mutate(
    the_accum_tsallis,
    diversity = exp_q(.data$entropy, q = q),
    .keep = "unused"
  )
  if (n_simulations > 0) {
    the_accum_hill <- dplyr::mutate(
      the_accum_hill,
      inf = exp_q(.data$inf, q = q),
      sup = exp_q(.data$sup, q = q)
    )
  }

  return(the_accum_hill)
}


#' @rdname accum_hill
#'
#' @export
accum_hill.abundances <- function(
    x,
    q = 0,
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

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  }
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  coverage_estimator <- match.arg(coverage_estimator)

  # Set levels if needed
  if (is.null(levels)) {
    sample_size <- max(
      rowSums(
        x[, !colnames(x) %in% non_species_columns]
      )
    )
    levels <- seq_len(sample_size)
  }
  # Apply accum_tsallis.numeric() to each site
  accum_hill_list <- apply(
    # Eliminate site and weight columns
    x[, !colnames(x) %in% non_species_columns],
    # Apply to each row
    MARGIN = 1,
    FUN = accum_hill.numeric,
    # Arguments
    q = q,
    levels = levels,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
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
  the_accum_hill <- tibble::tibble(
    site = rep(site_names, each = length(levels)),
    # Coerce the list returned by apply into a dataframe
    do.call(rbind.data.frame, accum_hill_list)
  )
  class(the_accum_hill) <- c("accumulation", class(the_accum_hill))

  return(the_accum_hill)
}
