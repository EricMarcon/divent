#' Shannon's Entropy of a Community
#' 
#' Estimate the entropy species from abundance or probability data.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#'
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities]
#' @param ... Unused.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#'
#' @return A tibble with the site names, the estimators used and the estimated entropy.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' ent_shannon(paracou_6_abd)
#' 
#' @name ent_shannon
NULL


#' @rdname ent_shannon
#'
#' @export
ent_shannon <- function(x, ...) {
  UseMethod("ent_shannon")
}


#' @rdname ent_shannon
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
ent_shannon.numeric <- function(
    x, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak",
                  "Grassberger2003", "Schurmann", "Bonachela", "Miller", "ZhangHz",
                  "naive"),
    level = NULL, 
    probability_estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    prob <- x[x > 0]
    ent_species <- -prob * log(prob)
    entropy <- sum(ent_species)
    if (as_numeric) {
      return(entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive", 
          entropy = entropy
        )
      )  
    }
  }
  
  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)
  # Probabilities
  prob <- abd / sample_size
  # Number of observed species
  s_obs <- length(abd)
  
  # Entropy of a vector of abundances ----
  if (is.null(level)) {
    ## Exit if x contains no or a single species ----
    if (s_obs < 2) {
      if (s_obs == 0) {
        if (as_numeric) {
          return(NA)
        } else {
          return(
            tibble::tibble_row(
              estimator = "No Species", 
              entropy = NA
            )
          )  
        }
      } else {
        if (as_numeric) {
          return(0)
        } else {
          return(
            tibble::tibble_row(
              estimator = "Single Species", 
              entropy = 0
            )
          )  
        }
      }
    } else {
      # Probabilities instead of abundances
      if (sample_size < 2) {
        warning("Entropy estimators can't apply to probability data. Estimator forced to 'naive'")
        estimator <- "naive"
      }
    }
    
    ## Naive estimator ----
    if (!is_integer_values(abd)) {
      warning("The estimator can't be applied to non-integer values.")
      estimator <- "naive"
    }
    if (estimator == "naive") {
      ent_species <- sum(ent_species)
      # Eliminate unobserved species. Useless because filtered before
      # ent_species[prob == 0] <- 0
      entropy <- sum(ent_species)
      if (as_numeric) {
        return(entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = entropy
          )
        )  
      }
    }
    
    ## Other estimators ----
    if (estimator == "Miller") {
      entropy <- ent_shannon.numeric(
          prob, 
          as_numeric = TRUE,
          check_arguments = FALSE
        ) +
        (s_obs - 1) / 2 / sample_size
      if (as_numeric) {
        return(entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = entropy
          )
        )  
      }
    }
    if (estimator == "ChaoShen" | estimator == "GenCov" | estimator == "Marcon") {
      sample_coverage <- coverage.numeric(abd, check_arguments = FALSE)$coverage
    }
    if (estimator == "ChaoShen") {
      prob_cov <- sample_coverage * prob
    }
    if (estimator == "GenCov" | estimator == "Marcon") {
      prob_cov <- probabilities.numeric(
        abd, 
        probability_estimator = probability_estimator, 
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha, 
        jack_max = jack_max,
        check_arguments = FALSE
      )
    } 
    if (estimator == "ChaoShen" | estimator == "GenCov" | estimator == "Marcon") {
      ent_ChaoShen <- -sum(prob_cov * log(prob_cov) / (1 - (1 - prob_cov)^sample_size))
    }
    if (estimator == "ChaoShen" | estimator == "GenCov") {
      return(tibble::tibble_row(estimator = estimator, entropy = ent_ChaoShen))
    } 
    if (estimator == "Grassberger" | estimator == "Marcon") {
      # (-1)^n is problematic for long vectors (returns NA for large values). 
      # It is replaced by 1 - n %%2 * 2 (abd is rounded if is not an integer)
      ent_Grassberger <- sum(
        abd / sample_size * 
        (log(sample_size) - digamma(abd) - (1 - round(abd) %%2 * 2) / (abd + 1))
      )
    }
    if (estimator == "Grassberger") {
      if (as_numeric) {
        return(ent_Grassberger)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = ent_Grassberger
          )
        ) 
      }
    }
    if (estimator == "Marcon") {
      entropy <- max(ent_ChaoShen, ent_Grassberger)
      if (as_numeric) {
        return(entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = entropy
          )
        ) 
      }
    }
    if (estimator == "Grassberger2003" | estimator == "Schurmann") {
      # Define a function to calculate the integral in the bias formula for each value of abd
      integral <- function(n, upper) {
        stats::integrate(function(t, n) t^(n - 1) / (1 + t), 0, upper, n) 
      }
    }
    if (estimator == "Grassberger2003") {
      integral_value <- unlist(
        vapply(
          abd, 
          integral, 
          FUN.VALUE = list(0.0, 0.0, 0, "", call("Integral", 0,0)), 
          upper = 1
        )["value",]
      )
    }
    if (estimator == "Schurmann") {
      integral_value <- unlist(
        vapply(
          abd, 
          integral, 
          FUN.VALUE=list(0.0, 0.0, 0, "", call("Integral", 0,0)), 
          upper = exp(-1/2)
        )["value",]
      )
    }
    if (estimator == "Grassberger2003" | estimator == "Schurmann") {
      entropy <- sum(
        abd / sample_size * 
        (digamma(sample_size) - digamma(abd) - (1 - abd %%2 * 2) * integral_value)
      )
      if (as_numeric) {
        return(entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = entropy
          )
        ) 
      }
    }
    if (estimator == "Holste" | estimator == "Bonachela") {
      seq_l <- seq_len(length(abd) + sample_size)
      inv_l <- 1 / seq_l
      cumul_l <- function(n) {sum(inv_l[n:length(inv_l)])}
      sum_inv_l <- vapply(seq_l, cumul_l, FUN.VALUE = 0) 
      if (estimator == "Holste") {
        ent_hb <- sum((abd + 1) / (length(abd) + sample_size) * sum_inv_l[abd + 2])
      } else {
        ent_hb <- sum((abd + 1) / (2 + sample_size) * sum_inv_l[abd + 2])
      }
      if (as_numeric) {
        return(ent_hb)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = ent_hb
          )
        ) 
      }
    }
    if (estimator == "ChaoJost") {
      # Calculate abundance distribution
      s_1 <- sum(abd == 1)
      s_2 <- sum(abd == 2)
      # Calculate A (Chao & Jost, 2015, eq. 6b)
      A <- chao_A(abd)
      # Chao, Wang & Jost 2013, eq. 7. Equals EntropyEstimation::Entropy.z(abd).
      ent_ChaoJost <- sum(
        abd / sample_size * (digamma(sample_size) - digamma(abd))
        )
      # Add Chao-Jost estimator to that of Zhang-Grabchak
      if (A != 1) {
        ent_part2 <- vapply(
          seq_len(sample_size-1), 
          function(r) 1 / r * (1 - A)^r, 
          FUN.VALUE = 0
        ) 
        ent_ChaoJost <- ent_ChaoJost + s_1 / sample_size * 
          (1 - A)^(1 -sample_size) * (-log(A) - sum(ent_part2))
      }
      if (as_numeric) {
        return(ent_ChaoJost)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = ent_ChaoJost
          )
        ) 
      }
    }
    if (estimator == "ZhangHz") {
      # Values of v in vector V
      V <- seq_len(sample_size - 1)
      # Weight part. Taken in log or goes to Inf for v > 1000
      # gamma cannot be used for large n, lgamma is preferred.
      lnw_v <- (V + 1) * log(sample_size) + lgamma(sample_size - V) - 
        lgamma(sample_size + 1) - log(V)
      # p_V_Ps is an array, containing (1 - p_s - j/n) for each species (lines) 
      # and all j from 0 to n-2. 
      # Because array indexation starts from 1 in R, j is replaced by j-1.
      p_V_prob <- outer(prob, V, function(p, j) {1 - p - (j - 1) / sample_size})
      # Useful values are products from j=0 to v, so prepare cumulative products
      p_V_prob <- t(apply(p_V_prob, 1, cumprod))
      # Sum of products, weighted by p_s
      sum_prod <- function(v) {sum(prob * p_V_prob[seq_along(prob), v])}
      # Apply sum_prod to all values of V. Use logs or w_v goes to Inf.
      entropy <- sum(exp(lnw_v + log(vapply(V, sum_prod, FUN.VALUE = 0))))
      if (as_numeric) {
        return(entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = entropy
          )
        ) 
      }
    }
    if (estimator == "UnveilC") {
      prob_unv <- probabilities.numeric(
        abd,
        estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = "Chao1", 
        as_numeric = TRUE,
        check_arguments = FALSE
      )
    }
    if (estimator == "UnveiliC") {
      prob_unv <- probabilities.numeric(
        abd,
        estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = "iChao1", 
        as_numeric = TRUE,
        check_arguments = FALSE
      )
    }
    if (estimator == "UnveilJ") {
      prob_unv <- probabilities.numeric(
        abd,
        estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = "jackknife", 
        jack_max = jack_max, 
        as_numeric = TRUE,
        check_arguments = FALSE
      )
    }
    if (estimator == "UnveilC" | estimator == "UnveiliC" | estimator == "UnveilJ") {
      # Naive estimator applied to unveiled distribution
      entropy <- -sum(prob_unv * log(prob_unv))
      if (as_numeric) {
        return(entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            entropy = entropy
          )
        )  
      }
    }
    
    # Not recognized
    warning("Correction was not recognized")
    return (NA)
  }

  # Entropy at a level ----
  # If Level is coverage, get size
  if (level < 1) {
    level <- coverage_to_size.numeric(
      abd, 
      sample_coverage = level,
      check_arguments = FALSE
    )
  }
  if (level == sample_size) {
    # No interpolation/extrapolation needed: return the observed entropy
    entropy <- -sum(prob * log(prob))
    if (as_numeric) {
      return(entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Sample", 
          entropy = entropy
        )
      )  
    }
  }
  if (level <= sample_size) {
    ## Interpolation ----
    abd_freq <- abd_freq_count(abd, level = level, check_arguments = FALSE)
    seq_l <- seq_len(level)/level
    entropy <- -(sum(seq_l * log(seq_l) * abd_freq$number_of_species))
    if (as_numeric) {
      return(entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Interpolation", 
          entropy = entropy
        )
      )  
    }
  } else {
    ## Extrapolation ----
    # Estimate the asymptotic entropy
    if (probability_estimator == "naive") {
      # Don't unveil the asymptotic distribution, use the asymptotic estimator
      ent_est <- ent_shannon(abd, estimator = estimator, check_arguments = FALSE)
    } else {
      # Unveil so that the estimation of H is similar to that of non-integer entropy
      prob_unv <- probabilities.numeric(
        abd,
        estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator, 
        jack_max = jack_max, 
        q = 1,
        as_numeric = TRUE,
        check_arguments = FALSE
      )
      ent_est <- -sum(prob_unv * log(prob_unv))
    }
    # Estimate observed entropy
    ent_obs <- -sum(prob * log(prob))
    # Interpolation
    entropy <- sample_size / level * ent_obs + 
      (level - sample_size) / level * ent_est
    if (as_numeric) {
      return(entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = probability_estimator, 
          entropy = entropy
        )
      )  
    }
  }
}


#' @rdname ent_shannon
#'
#' @param gamma If `TRUE`, $\gamma$ diversity, i.e. diversity of the metacommunity, is computed.
#' 
#' @export
ent_shannon.abundances <- function(
    x,
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak"),
    level = NULL, 
    probability_estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  if (gamma) {
    return(
      ent_shannon.numeric(
        metacommunity(x, as_numeric = TRUE, check_arguments = FALSE),
        # Arguments
        estimator = estimator,
        level = level, 
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha, 
        jack_max = jack_max, 
        as_numeric = FALSE,
        check_arguments = FALSE
      )
    )
  } else {
    # Apply ent_shannon.numeric() to each site
    ent_shannon_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% c("site", "weight"))], 
      # Apply to each row
      MARGIN = 1,
      FUN = ent_shannon.numeric,
      # Arguments
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max, 
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    return(
      # Make a tibble with site, estimator and richness
      tibble::tibble(
        # Do not assume column site exists
        x[colnames(x) == "site"],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, ent_shannon_list)
      )
    )
  }
}