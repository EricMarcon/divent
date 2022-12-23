#' Tsallis Entropy of a Community
#' 
#' Estimate the entropy of species from abundance or probability data.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#' See [div_hill] for estimators.
#' 
#' Entropy can be estimated at a specified level of interpolation or 
#' extrapolation, either a chosen sample size or sample coverage 
#' \insertCite{Chao2014}{divent}, rather than its asymptotic value.
#' See [ent_accum] for details.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @return A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' ent_tsallis(paracou_6_abd, q = 2)
#' 
#' # At 80% coverage
#' ent_tsallis(paracou_6_abd, level=0.8)
#' 
#' @name ent_tsallis
NULL


#' @rdname ent_tsallis
#'
#' @export
ent_tsallis <- function(x, q = 1, ...) {
  UseMethod("ent_tsallis")
}


#' @rdname ent_tsallis
#'
#' @param estimator An estimator of entropy. 
#' 
#' @export
ent_tsallis.numeric <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    sample_coverage = NULL,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    if (!is.null(level)) stop("Entropy can't be estimated at a level without the abundance of species.")
    # Probabilities sum to 1, allowing rounding error
    prob <- x[x > 0]
    # Avoid big numbers with this rather than prob * ln_q(1/prob, q)
    ent_species <- -prob^q * ln_q(prob, q = q)
    the_entropy <- sum(ent_species)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive", 
          order = q,
          entropy = the_entropy
        )
      )  
    }
  }
  
  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)
  # Number of observed species
  s_obs <- length(abd)
  
  # Entropy of a vector of abundances ----
  if (is.null(level)) {
    ## Exit if x contains no or a single species ----
    if (length(abd) < 2) {
      if (length(abd) == 0) {
        if (as_numeric) {
          return(NA)
        } else {
          return(
            tibble::tibble_row(
              estimator = "No Species", 
              order = q,
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
              order = q,
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
    
    ## Metacommunity estimation ----
    # abd may not be integers, then sample_coverage is given
    # sample_coverage is between 0 and 1 (by check_arguments), sum(abd) must be an integer.
    # estimator may be ChaoShen or Marcon (max(ChaoShen, Grassberger))
    if (
        !is.null(sample_coverage) & 
        is_integer_values(sample_size) & 
        (estimator == "ChaoShen" | estimator == "Marcon")) {
      
      cp <- sample_coverage * abd / sample_size
      chao_shen <- -sum(cp^q * ln_q(cp, q = q) /(1 - (1 - cp)^sample_size))
      if (estimator == "Marcon") {
        # Calculate Grassberger's correction
        if (q == 1) {
          grassberger <- sum(
            abd / sample_size * (log(sample_size) - digamma(abd) - 
            (1 - round(abd) %% 2 * 2) / (abd + 1))
          )
        } else {
          grassberger <- (1 - sample_size^(-q) * sum(e_n_q(abd, q = q))) / (q - 1)
        }
      } else grassberger <- 0
      # Take the max
      if (chao_shen > grassberger) {
        if (as_numeric) {
          return(chao_shen)
        } else {
          return(
            tibble::tibble_row(
              estimator = "ChaoShen", 
              order = q,
              entropy = chao_shen
            )
          )  
        }
      } else {
        if (as_numeric) {
          return(grassberger)
        } else {
          return(
            tibble::tibble_row(
              estimator = "Grassberger", 
              order = q,
              entropy = grassberger
            )
          )  
        }
      }
    }
    
    ## Naive estimator ----
    if (!is_integer_values(abd)) {
      warning("The estimator can't be applied to non-integer values.")
      estimator <- "naive"
    }
    if (estimator == "naive") {
      prob <- abd / sample_size
      # Avoid big numbers with this rather than prob * ln_q(1/prob, q)
      # Eliminate unobserved species. Useless because filtered before
      # ent_species <- -prob^q * ln_q(prob, q)
      # ent_species[prob == 0] <- 0
      the_entropy <- -sum(prob^q * ln_q(prob, q))
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = the_entropy
          )
        ) 
      }
    }
    
    # Common code for ZhangGrabchak. Useless if EntropyEstimation is used.
    # if (estimator == "ZhangGrabchak" | estimator == "ChaoWangJost" | estimator == "ChaoJost") {
    #   prob <- abd / sample_size
    #   V <- seq_len(sample_size - 1)
    #   # p_V_abd is an array, containing (1 - (s_obs - 1) / (sample_size - j)) 
    #   # for each species (lines) and all j from 1 to sample_size - 1
    #   p_V_abd <- outer(abd, V, function(abd, j) 1- (abd - 1) / (sample_size - j))
    #   # Useful values are products from j = 1 to v, so prepare cumulative products
    #   p_V_abd <- t(apply(p_V_abd, 1, cumprod))
    #   # Sum of products weighted by w_v
    #   sum_prod <- function(s) {
    #     used_v <- seq_len(sample_size - abd[s])
    #     return(sum(w_v[used_v] * p_V_abd[s, used_v]))
    #   }
    # }
    # Specific code for q=1
    # # Weights
    # w_v <- 1/V
    # the_entropy <- sum(prob * vapply(seq_along(abd), S_v, 0))
    
    
    
    ## Shannon ----
    if (q == 1) {
      return(
        ent_shannon.numeric(
          abd, 
          estimator = estimator, 
          probability_estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator,
          jack_alpha  = jack_alpha, 
          jack_max = jack_max,
          as_numeric = as_numeric,
          check_arguments = FALSE
        )
      )
    }
    
    ## Not Shannon ----
    if (estimator == "ZhangGrabchak" | estimator == "ChaoJost") {
      # Weights. Useless here if EntropyEstimation is used, but weights are necessary for ChaoJost
      # i <- seq_len(abd)
      # w_vi <- (i - q) / i
      # w_v <- cumprod(w_vi)
      # ZhangGrabchak <- sum(prob * vapply(seq_along(abd), sum_prod, 0)) / (1 - q)
      # Use EntropyEstimation instead
      if (q==0) {
        ent_ZhangGrabchak <- s_obs-1 
      } else {
        ent_ZhangGrabchak <- EntropyEstimation::Tsallis.z(abd, q)
      }
      # ZhangGrabchak stops here, but ChaoWangJost adds an estimator of the bias
      if (estimator == "ZhangGrabchak") {
        if (as_numeric) {
          return(ent_ZhangGrabchak)
        } else {
          return(
            tibble::tibble_row(
              estimator = estimator, 
              order = q,
              entropy = ent_ZhangGrabchak
            )
          )  
        }
      }
      # Calculate abundance distribution
      s_1 <- sum(abd == 1)
      s_2 <- sum(abd == 2)
      # Calculate A (Chao & Jost, 2015, eq. 6b)
      A <- chao_A(abd)
      # Eq 7d in Chao & Jost (2015). 
      # Terms for r in seq_len(sample_size-1) equal (-1)^r * w_v[r] * (A-1)^r. 
      # w_v is already available from ZhangGrabchak
      i <- seq_len(sample_size)
      # Weights: here only if EntropyEstimation is used. Else, they have been calculated above.
      w_vi <- (i - q)/i
      w_v <- cumprod(w_vi)
      if (A == 1) {
        # The general formula of Eq 7d has a 0/0 part that must be forced to 0
        bias_chao_jost <- 0
      } else {
        eq7d_sum <- vapply(
          seq_len(sample_size - 1), 
          function(r) {
            w_v[r] * (1 - A)^r
          }, 
          FUN.VALUE = 0
        )
        # Calculate the estimator of the bias. 
        # eq7d_sum contains all terms of the sum except for r=0: the missing term equals 1.
        # The bias in Chao & Jost (2015) is that of the Hill number. 
        # It must be divided by 1-q to be applied to entropy.
        bias_chao_jost <- (
          s_1 / sample_size * (1 - A)^(1 - sample_size) * (A^(q - 1) - sum(eq7d_sum) - 1)
          ) / (1 - q)
      }
      the_entropy <- ent_ZhangGrabchak + bias_chao_jost
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = the_entropy
          )
        )  
      }
    }
    if (estimator == "ChaoShen" | estimator == "Marcon") {
      sample_coverage <- coverage.numeric(
        abd, 
        as_numeric = TRUE,
        check_arguments = FALSE
      )
      prob_cov <- sample_coverage * abd / sample_size
    }
    if (estimator == "GenCov") {
      if (probability_estimator == "naive") {
        warning("'probability_estimator' can't be 'naive' with estimator 'GenCov'. It is forced to 'Chao2015'.")
      }
      if (unveiling != "none") {
        warning("'unveiling' is forced to 'none' with estimator 'GenCov'.")
      }
      prob_cov <- probabilities.numeric(
        abd, 
        estimator = probability_estimator, 
        unveiling = "none",
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha, 
        jack_max = jack_max,
        as_numeric = TRUE,
        check_arguments = FALSE
      )
    } 
    if (estimator == "ChaoShen" | estimator == "Marcon" | estimator == "GenCov") {
      ent_cov <- -sum(prob_cov^q * ln_q(prob_cov, q) / (1 - (1 - prob_cov)^sample_size))
    }
    if (estimator == "ChaoShen" | estimator == "GenCov") {
      if (as_numeric) {
        return(ent_cov)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = ent_cov
          )
        )  
      }
    }
    if (estimator == "Grassberger" | estimator == "Marcon") {
      ent_Grassberger <- (1 - sample_size^(-q) * sum(e_n_q(abd, q))) / (q - 1)
    }
    if (estimator == "Grassberger") {
      if (as_numeric) {
        return(ent_Grassberger)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = ent_Grassberger
          )
        )  
      }
    }
    if (estimator == "Marcon") {
      the_entropy <- max(ent_cov, ent_Grassberger)
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = the_entropy
          )
        )  
      }
    }
    if (estimator == "Holste") {
      the_entropy <- 1 / (1 - q) * 
        (beta(s_obs + sample_size, q) * sum(1 / beta(abd + 1, q)) - 1)
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = the_entropy
          )
        )  
      }
    } 
    if (estimator == "Bonachela") {
      the_entropy <- 1 / (1 - q) *
        (beta(2 + sample_size, q) * sum(1 / beta(abd + 1, q)) - 1)
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = the_entropy
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
      the_entropy <- -sum(prob_unv^q * ln_q(prob_unv, q))
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = q,
            entropy = the_entropy
          )
        )
      }
    }
    
    warning("estimator was not recognized")
    return (NA)
  }
  
  # Entropy at a level ----

  # Single species. 
  if (s_obs == 1 | level == 1) {
    if (as_numeric) {
      return(0)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Single Species", 
          order = q,
          level = 1,
          entropy = 0
        )
      )  
    }
  }
  
  # If level is coverage, get size
  if (level < 1) {
    level <- coverage_to_size.numeric(
      abd, 
      sample_coverage = level,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
  }
  if (level == sample_size) {
    # No interpolation/extrapolation needed: return the observed entropy
    the_entropy <- -sum(prob^q * ln_q(prob, q))
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Sample", 
          order = q,
          level = level,
          entropy = the_entropy
        )
      )  
    }
  }
  
  ## Integer q ----
  if (q==0) {
    # Richness - 1. Same result as general formula but faster
    the_richness <- div_richness.numeric(
      abd, 
      # Unused
      estimator = richness_estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max, 
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    # Calculate entropy and return
    if (as_numeric) {
      return(the_richness$diversity - 1)
    } else {
      the_richness <- dplyr::mutate(
        the_richness,
        entropy = .data$diversity - 1,
        .keep = "unused")
      return(the_richness)  
    }
  } 
  if (q==1) {
    # Shannon. General formula is not defined at q=1
    return(
      ent_shannon(
        x, 
        estimator = estimator,
        level = level, 
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha, 
        jack_max = jack_max,
        as_numeric  = as_numeric,
        check_arguments = FALSE
      )
    )
  } 
  if (q==2) {
    # Simpson.
    return(
      ent_simpson(
        x, 
        # No estimator needed. Avoid raising an error with an unknown one.
        estimator = "naive",
        level = level, 
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha, 
        jack_max = jack_max,
        as_numeric  = as_numeric,
        check_arguments = FALSE
      )
    )
  }
  
  if (level <= sample_size) {
    ## Interpolation ----
    # Obtain Abundance Frequence Count
    abd_freq <- abd_freq_count(abd, level = level, check_arguments = FALSE)
    # Calculate entropy (Chao et al., 2014, eq. 6)
    the_entropy <- (sum(((seq_len(level)) / level)^q * abd_freq$number_of_species) - 1)/(1 - q)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Interpolation", 
          order = q,
          level = level,
          entropy = the_entropy
        )
      )  
    }
  } else {
    ## Extrapolation ----
    # Unveil the full distribution that rarefies to the observed entropy
    prob_unv <- probabilities.numeric(
      abd,
      estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator, 
      jack_max = jack_max, 
      q = q,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    # Obtain Abundance Frequence Count at level (Chao et al., 2014, eq. 5)
    s_level <- vapply(
      seq_len(level), 
      function(nu) {
        sum(
          exp(
            lchoose(level, nu) + nu * log(prob_unv) + (level - nu) * log(1 - prob_unv)
          )
        )
      }, 
      FUN.VALUE=0
    )
    # Estimate entropy (Chao et al., 2014, eq. 6)
    the_entropy <- (sum((seq_len(level) / level)^q * s_level) - 1) / (1 - q)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = richness_estimator, 
          order = q,
          level = level,
          entropy = the_entropy
        )
      )  
    }
  }
}


#' @rdname ent_tsallis
#'
#' @export
ent_tsallis.species_distribution <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  if (gamma) {
    return(
      ent_gamma(
        x = x,
        q = q,
        estimator = estimator,
        level = level, 
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha,
        jack_max = jack_max
      )
    )
  } else {
    # Apply ent_tsallis.numeric() to each site
    ent_tsallis_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% non_species_columns)], 
      # Apply to each row
      MARGIN = 1,
      FUN = ent_tsallis.numeric,
      # Arguments
      q = q,
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
      # Make a tibble with site, estimator and entropy
      tibble::tibble(
        # Restore non-species columns
        x[colnames(x) %in% non_species_columns],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, ent_tsallis_list)
      )
    )
  }
}


#' Gamma entropy of a metacommunity
#' 
#' Build the metacommunity and check that abundances are integers.
#' If they are not (due to weights of communities) then use a fallback estimator:
#' "ChaoShen" requires the sample coverage of the assemblage of sites.
#' "Grassberger" accepts non integer abundances.
#' "Marcon" combines both.
#' 
#' `ent_tsallis.numeric` contains the only implementation of this estimation.
#' i.e., `ent_shannon` can't be used but `ent_tsallis.numeric`with `q=1` will work fine.
#'
#' @return A tibble with the estimator used and the estimated entropy.
#' @noRd
ent_gamma <- function(
    x,
    q,
    estimator,
    level,
    probability_estimator,
    unveiling,
    richness_estimator,
    jack_alpha,
    jack_max) {
  
  # Build the metacommunity
  abd <- metacommunity(
    x, 
    as_numeric = TRUE, 
    check_arguments = FALSE
  )
  if (is_integer_values(abd)) {
    # Sample coverage is useless
    sample_coverage <- NULL
  } else {
    # Non integer values in the metacommunity. 
    # Calculate the sample coverage and change the estimator.
    sample_coverage <- coverage.numeric(
      colSums(x[, !(colnames(x) %in% non_species_columns)]),
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    if (!estimator %in% c("Marcon", "ChaoShen")) {
      estimator <- "Marcon"
    }
  }
  return(
    ent_tsallis.numeric(
      abd,
      # Arguments
      q = q,
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max,
      sample_coverage = sample_coverage,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
  )
}