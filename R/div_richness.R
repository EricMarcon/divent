#' Number of Species of a Community
#' 
#' Estimate the number of species from abundance or probability data.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#' Chao's correction techniques are from \insertCite{Chao2014;textual}{divent} 
#' and \insertCite{Chiu2014a;textual}{divent}. 
#' The Jackknife estimator is calculated by a straight adaptation of the code 
#' by Ji-Ping Wang (jackknife in package **SPECIES**). 
#' The optimal order is selected according to 
#' \insertCite{Burnham1978,Burnham1979;textual}{divent}. 
#' Many other estimators are available elsewhere, the ones implemented here are 
#' necessary for other entropy estimations.
#' 
#' Richness can be estimated at a specified `level` of interpolation or 
#' extrapolation, either a chosen sample size or sample coverage 
#' \insertCite{Chiu2014a}{divent}, rather than its asymptotic value. 
#' Extrapolation relies on the estimation of the asymptotic richness. 
#' If `probability_estimator` is "naive", then the asymptotic estimation of 
#' richness is made using the chosen `estimator`, else the asymptotic 
#' distribution of the community is derived and its estimated richness adjusted 
#' so that the richness of a sample of this distribution of the size of the 
#' actual sample has the richness of the actual sample.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities]
#' @param ... Unused.
#' @param gamma If `TRUE`, gamma diversity is calculated.
#' The metacommunity if built by combining the community abundances with respect to their weight.
#'
#' @return A tibble with the site names, the estimators used and the estimated numbers of species.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' div_richness(paracou_6_abd)
#' 
#' @name div_richness
NULL


#' @rdname div_richness
#'
#' @export
div_richness <- function(x, ...) {
  UseMethod("div_richness")
}


#' @rdname div_richness
#'
#' @param estimator An estimator of richness to evaluate the total number of species. 
#' @param level The level of interpolation or extrapolation. 
#' It may be a sample size (an integer) or a sample coverage 
#' (a number between 0 and 1).
#' The asymptotic `estimator` is used in extrapolation 
#' (i.e. a `level` greater than the sample size).
#' @param probability_estimator A string containing one of the possible estimators
#' of the probability distribution (see [probabilities]). 
#' Used only by the estimator of richness "rarefy".
#' @param unveiling A string containing one of the possible unveiling methods 
#' to estimate the probabilities of the unobserved species (see [probabilities]).
#' Used only by the estimator of richness "rarefy".
#' 
#' @export
div_richness.numeric <- function(
    x, 
    estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10, 
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    as_numeric = FALSE,
    ..., 
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Diversity of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    if (!is.null(level)) stop("Richness can't be estimated at a level without the abundance of species.")
    # Probabilities sum to 1, allowing rounding error
    return(
      tibble::tibble_row(
        estimator = "naive", 
        order = 0,
        diversity = sum(x > 0)
      )
    )
  }
  
  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)
  # Number of observed species
  s_obs <- length(abd)

  # Diversity of a vector of abundances ----
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
              order = 0,
              diversity = NA
            )
          )
        }
      } else {
        if (as_numeric) {
          return(1)
        } else {
          return(
            tibble::tibble_row(
              estimator = "Single Species", 
              order = 0,
              diversity = 1
            )
          )
        }
      }
    } else {
      # Probabilities instead of abundances
      if (sample_size < 2) {
        warning("Richness estimators can't apply to probability data. Estimator forced to 'naive'")
        estimator <- "naive"
      }
    }
    
    ## Naive estimator ----
    if (estimator == "naive") {
      if (as_numeric) {
        return(s_obs)
      } else {
        return(
          tibble::tibble_row(
            estimator = "naive", 
            order = 0,
            diversity = s_obs
          )
        )
      }
    }
    
    ## Rarefaction estimator ----
    if (estimator == "rarefy") {
      the_richness <- length(
        probabilities.numeric(
          abd, 
          estimator = probability_estimator, 
          unveiling = unveiling, 
          richness_estimator = "rarefy", 
          q = 0, 
          check_arguments = FALSE
        )
      )
      if (as_numeric) {
        return(the_richness)
      } else {
        return(
          tibble::tibble_row(
            estimator = "rarefy", 
            order = 0,
            diversity = the_richness
          )
        )
      }
    }
    
    ## Abundance frequency count ----
    abd_freq <- abd_freq_count(abd)
    s_1 <- as.integer(abd_freq[abd_freq[, 1] == 1, 2])
    s_2 <- as.integer(abd_freq[abd_freq[, 1] == 2, 2])
    
    
    ## Chao1 and iChao1 ----
    if ((estimator == "Chao1") | (estimator == "iChao1")) {
      if (is.na(s_1)) {
        s_1 <- 0
      }
      if (is.na(s_2)) {
        s_0 <- (sample_size - 1) / sample_size * s_1 * (s_1 - 1) / 2
        s_2 <- 0
      } else {
        s_0 <- (sample_size - 1) / sample_size * s_1 * s_1 / 2 / s_2
      }
    }
    ### Chao1 ----
    if (estimator == "Chao1") {
      the_richness <- s_obs + s_0
      if (as_numeric) {
        return(the_richness)
      } else {
        return(
          tibble::tibble_row(
            estimator = "Chao1", 
            order = 0,
            diversity = the_richness
          )
        )
      }
    }
    ### iChao1 ----
    if (estimator == "iChao1") {
      s_3 <- as.integer(abd_freq[abd_freq[, 1] == 3, 2])
      s_4 <- as.integer(abd_freq[abd_freq[, 1] == 4, 2])
      if (is.na(s_3)) {
        s_3 <- 0
      }
      if (is.na(s_4)) {
        s_4 <- 1
      }
      s_0_iChao <- s_3 / 4 / s_4 * max(s_1 - s_2 * s_3 / 2 / s_4, 0)
      the_richness <- s_obs + s_0 + s_0_iChao
      if (as_numeric) {
        return(the_richness)
      } else {
        return(
          tibble::tibble_row(
            estimator = "iChao1", 
            order = 0,
            diversity = the_richness
          )
        )
      }
    }
    
    ## jackknife ----
    if (estimator == "jackknife") {
      # Adapted from jackknife in SPECIES
      # Max possible order
      k <- min(nrow(abd_freq) - 1, jack_max)
      if (k == 0) {
        # No optimisation possible. Return "jackknife of order 0".
        k_smallest <- 1
        s_jack <- s_obs
      } else {
        # Max number of individuals
        m <- max(abd_freq[, 1])
        # Complete the abundance frequency count for all counts between 1 and m
        n_temp <- cbind(seq_len(m), rep(0, m))
        n_temp[abd_freq$abundance , 2] <- abd_freq$number_of_species
        abd_freq <- n_temp
        # Prepare a matrix with k+1 rows and 5 columns
        gene <- matrix(0, nrow = k + 1, ncol = 5)
        gene[1, 1] <- s_obs
        for (i in seq_len(k)) {
          gene[i + 1, 1] <- s_obs
          gene[i + 1, 4] <- s_obs
          for (j in seq_len(i)) {
            gene[i + 1, 1] <- gene[i + 1, 1] + (-1)^(j + 1) * 
              2^i * stats::dbinom(j, i, 0.5) * abd_freq[j, 2]
            gene[i + 1, 4] <- gene[i + 1, 4] + (-1)^(j + 1) * 
              2^i * stats::dbinom(j, i, 0.5) * abd_freq[j, 2] * prod(seq_len(j))
          }
          gene[i + 1, 2] <- -gene[i + 1, 1]
          for (j in seq_len(i)) {
            gene[i + 1, 2] <- gene[i + 1, 2] + ((-1)^(j + 1) *
              2^i * stats::dbinom(j, i, 0.5) + 1)^2 * abd_freq[j, 2]
          }
          gene[i + 1, 2] <- gene[i + 1, 2] + sum(abd_freq[(i + 1):nrow(abd_freq), 2])
          gene[i + 1, 2] <- sqrt(gene[i + 1, 2])
        }
        if (k > 1) {
          for (i in 2:k) {
            gene[i, 3] <- -(gene[i + 1, 1] - gene[i, 1])^2/(s_obs - 1)
            for (j in seq_len(i - 1)) {
              gene[i, 3] <- gene[i, 3] + 
                (
                  (-1)^(j + 1) * 2^(i) * stats::dbinom(j, i, 0.5) - 
                  (-1)^(j + 1) * 2^(i - 1) * stats::dbinom(j, i - 1, 0.5)
                )^2 * 
                abd_freq[j, 2] * s_obs/(s_obs - 1)
            }
            gene[i, 3] <- gene[i, 3] + abd_freq[i, 2] * s_obs/(s_obs - 1)
            gene[i, 3] <- sqrt(gene[i, 3])
            gene[i, 5] <- (gene[i + 1, 1] - gene[i, 1])/gene[i, 3]
          }
        }
        # Threshold for Burnham and Overton's test
        coe <- stats::qnorm(1 - jack_alpha/2, 0, 1)
        # Which orders pass the test?
        orders <- (gene[2:(k + 1), 5] < coe)
        if (sum(orders, na.rm=TRUE) == 0) {
          # If none, keep the max value of k (+1 because jack1 is in line 2)
          k_smallest <- k + 1
          s_jack <- gene[k_smallest, 1]
          # Estimated standard error
          sej <- gene[k_smallest, 2]
        }
        else {
          # Else, keep the smallest value (+1 because jack1 is in line 2)
          k_smallest <- which(orders)[1] + 1
          # Estimated value
          s_jack <- gene[k_smallest, 1]
          # Estimated standard error
          sej <- gene[k_smallest, 2]
        }
      }
      if (as_numeric) {
        return(s_jack)
      } else {
        return(
          tibble::tibble_row(
            estimator = paste("Jackknife", k_smallest - 1),
            order = 0,
            diversity = s_jack
          )
        )
      }
    }
  }

  # Diversity at a level ----
  # If level is coverage, get size
  if (level < 1) {
    level <- coverage_to_size.numeric(
      abd, 
      sample_coverage = level,
      check_arguments = FALSE
    )
  }
  if (level == sample_size) {
    # No interpolation/extrapolation needed: return the observed number of species
    if (as_numeric) {
      return(s_obs)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Sample", 
          order = 0,
          level = level,
          diversity = s_obs
        )
      )
    }
  }
  if (level <= sample_size) {
    ## Interpolation ----
    the_richness <- s_obs - 
      sum(
        exp(lchoose(sample_size - abd, level) - lchoose(sample_size, level))
      )
    if (as_numeric) {
      return(the_richness)
    } else {
      return(
        tibble::tibble_row(
          estimator = "SAC",
          order = 0,
          level = level,
          diversity = the_richness
        )
      )  
    }
  } else {
    ## Extrapolation ----
    s_1 <-  sum(x == 1)
    if (s_1) {
      # Estimate the number of unobserved species
      if (probability_estimator == "naive") {
        # Don't unveil the asymptotic distribution, use the asymptotic estimator
        s_0 <- div_richness.numeric(
          abd, 
          estimator=estimator, 
          jack_alpha = jack_alpha, 
          jack_max = jack_max, 
          as_numeric = TRUE,
          check_arguments = FALSE
        ) - s_obs
      } else {
        # Unveil so that the estimation of richness is similar to that of non-integer entropy
        prob_s_0 <- probabilities.numeric(
          abd, 
          estimator = probability_estimator, 
          unveiling = unveiling, 
          richness_estimator = estimator, 
          q = 0,
          as_numeric = TRUE,
          check_arguments = FALSE
        )
        s_0 <- length(prob_s_0) - s_obs
      }
      the_richness <- s_obs + s_0 * (
        1 - (1 - s_1 / (sample_size * s_0 + s_1))^(level - sample_size)
      )
    } else {
      # No singleton
      the_richness <- s_obs
    }
    if (as_numeric) {
      return(the_richness)
    } else {
      return(
        tibble::tibble_row(
          estimator = estimator, 
          order = 0,
          level = level,
          diversity = the_richness
        )
      )  
    }
  }
}


#' @rdname div_richness
#'
#' @param gamma If `TRUE`, \eqn{\gamma} diversity, i.e. diversity of the metacommunity, is computed.
#' 
#' @export
div_richness.species_distribution <- function(
    x, 
    estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10, 
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    gamma = FALSE,
    ..., 
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  if (gamma) {
    # Choose the estimator of entropy
    if (estimator == "rarefy") stop ("The 'rarefy' estimator is not implemented for gamma diversity.")
    # Entropy estimators rely on richness estimators
    ent_estimator <- switch(
      estimator,
      jackknife = "UnveilJ", 
      iChao1 = "UnveiliC", 
      Chao1 = "UnveilC", 
      naive = "naive"
    )
    # Calculate gamma entropy of order 0
    ent_0 <- ent_gamma(
      x = x,
      q = 0,
      estimator = ent_estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = estimator,
      jack_alpha  = jack_alpha,
      jack_max = jack_max
    )
    # Calculate diversity
    ent_0 <- dplyr::mutate(
      ent_0, 
      diversity = .data$entropy + 1, 
      .keep = "unused")
    # return the richness
    return(ent_0)
  } else {
    # Apply div_richness.numeric() to each site
    div_richness_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% c("site", "weight"))], 
      # Apply to each row
      MARGIN = 1,
      FUN = div_richness.numeric,
      # Arguments
      estimator = estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max, 
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    return(
      # Make a tibble with site, estimator and richness
      tibble::tibble(
        # Do not assume column site exists
        x[colnames(x) == "site"],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, div_richness_list)
      )
    ) 
  }
}
