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
#' # TODO
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
#' @param estimator An estimator of entropy. 
#' @param level The level of interpolation or extrapolation. 
#' It may be a chosen sample size (an integer) or a sample coverage 
#' (a number between 0 and 1). 
#' Richness extrapolation require its asymptotic estimation depending on the 
#' choice of the `estimator`.
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
#' @param as_numeric If `TRUE`, a number is returned rather than a tibble.
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
    as_numeric  = FALSE,
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
  
  # Level == Sample Size ----
  if (any(levels == sample_size)) {
    i <- which(levels == sample_size)
    ent_level[i] <- ent_tsallis.numeric(
      prob, 
      q = q,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
  }
  
  # Extrapolation ----
  # Don't use Tsallis for speed: probabily extrapolation should be run once only.
  levels_extrap <- levels[levels > sample_size]
  prob_unv <- NULL
  if (length(levels_extrap)) {
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
      if(show_progress & interactive()) 
        utils::setTxtProgressBar(pgb, length(levels))
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
          if(show_progress & interactive()) utils::setTxtProgressBar(pgb, length(levels))
        } else {
          # General case: q is not 0, 1 or 2 ----
          for(level in levels_extrap) {
            # Abundance frequence count at Level (Chao et al., 2014, eq. 5)
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
            i <- which(levels == Level)
            ent_level[i] <- (sum((seq_len(level) / level)^q * s_nu) - 1) / (1 - q)
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
    ent_sim <- matrix(0, nrow = length(levels), ncol = 2)
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
    for(Level in levels) {
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
        checkArguments=FALSE
      )
      i <- which(levels == level)
      # Store quantiles
      ent_sim[i, ] <- stats::quantile(ent_sim, probs = c(alpha / 2, 1 - alpha / 2))
      if (show_progress & interactive()) utils::setTxtProgressBar(pgb, i)
    }
  }
  close(pgb)

  # Format the result ----
  if (n_simulations) {
    accumulation <- tibble::tibble(
      level = levels,
      entropy = ent_level,
      inf = ent_sim[, 1],
      sup = ent_sim[, 2]
    )
  } else {
    accumulation <- tibble::tibble(
      level = levels,
      entropy = ent_level
    )
  }

  # Return actual values as attributes
  attr(accumulation, "sample_size") <- sample_size
  if (any(levels == sample_size)) {
    attr(accumulation, "entropy") <- ent_level[which(levels == sample_size)]
  } else {
    attr(accumulation, "entropy") <- ent_tsallis.numeric(
      prob, 
      q=q,
      check_arguments = FALSE
    )
  }
  class(accumulation) <- c("ent_accumulation", "accumulation", class(accumulation))
  
  return(accumulation)
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
      colSums(
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
    
  return(
    # Make a tibble with site, estimator and richness
    tibble::tibble(
      site = rep(site_names, each = length(levels)),
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, ent_accum_list)
    )
  )
}


#' @rdname div_accum
#'
#' @export
div_accum <- function(x, ...) {
  UseMethod("div_accum")
}
