#' Similarity-Based Entropy of a Community
#' 
#' Estimate the entropy of species from abundance or probability data and a
#' similarity matrix between species.
#' Several estimators are available to deal with incomplete sampling. 
#' Bias correction requires the number of individuals.
#' 
#' All species of the `species_distribution` must be found in the matrix of 
#' `similarities` if it is named.
#' If it is not, its size must equal the number of species.
#' Then, the order of species is assumed to be the same as that of the
#' `species_distribution`.
#' 
#' Similarity-Based entropy can't be interpolated of extrapolated as of the
#' state of the art.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' If it is a numeric vector, then its length must equal the dimensions of the
#' `similarities` matrix: species are assumed to be in the same order.
#' @param ... Unused.
#'
#' @return A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Entropy of each community
#' ent_similarity(paracou_6_abd, q = 2)
#' # gamma entropy
#' ent_similarity(paracou_6_abd, q = 2, gamma = TRUE)
#' 
#' @name ent_similarity
NULL


#' @rdname ent_similarity
#'
#' @export
ent_similarity <- function(x, similarities, q = 1, ...) {
  UseMethod("ent_similarity")
}


#' @rdname ent_similarity
#'
#' @param estimator An estimator of entropy. 
#' 
#' @export
ent_similarity.numeric <- function(
    x, 
    similarities = diag(length(x)),
    q = 1, 
    estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang", 
                  "UnveilC", "UnveiliC", "naive"),
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
     jack_alpha  = 0.05, 
    jack_max = 10,
    sample_coverage = NULL,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Check that dimensions correspond
  if (length(x) != ncol(similarities)) {
    stop("The length of 'x' must be equal to the dimension of the similarities.")
  }

  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error.
    if (q == 1) {
      # Limit value
      the_entropy <- as.numeric(-x %*% log(ordinariness))
    } else {
      the_entropy <- as.numeric((1 - (x %*% (ordinariness^(q - 1)))) / (q - 1))
    }
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
  similarities <- similarities[x > 0, x > 0]
  # Sample size
  sample_size <- sum(abd)
  # Calculate ordinariness
  prob <- abd / sample_size
  ordinariness <- as.numeric(similarities %*% prob)
  # Number of observed species
  s_obs <- length(abd)
  
  
  # Entropy of a vector of abundances ----
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
  if (!is.null(sample_coverage)) {
    # Force estimator to ChaoShen
    estimator <- "ChaoShen"
  } else {
    # Calculate sample coverage
    sample_coverage <- coverage.numeric(
      abd, 
      estimator = coverage_estimator,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
  }
  
  ## Naive estimator ----
  if (!is_integer_values(abd)) {
    warning("The estimator can't be applied to non-integer values.")
    estimator <- "naive"
  }
  if (estimator == "naive") {
    if (q == 1) {
      # Limit value
      the_entropy <- as.numeric(-prob %*% log(ordinariness))
    } else {
      the_entropy <- as.numeric((1 - (prob %*% (ordinariness^(q - 1)))) / (q - 1))
    }
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
  
  ## Unveiled estimator ----
  if (estimator %in% c("UnveilJ", "UnveilC", "UnveiliC")) {
    prob_unv <- probabilities.numeric(
      abd,
      estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = switch(
        estimator,
        UnveilJ = "jackknife", 
        UnveilC = "Chao1", 
        UnveiliC = "iChao1"
      ), 
      jack_alpha = jack_alpha,
      jack_max = jack_max,
      coverage_estimator = coverage_estimator,
      q = q,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    s_est <- length(prob_unv)
  }
  
  # Calculate the average similarity
  Z <- similarities
  diag(Z) <- NA
  sim_mean <- mean(Z, na.rm = TRUE)
  # Tune probabilities
  prob_cov <- prob * sample_coverage
  # Tune the ordinariness of observed species, i.e.
  # similarity with observed species is multiplied by coverage
  # and unobserved species add (1 - coverage) * average similarity 
  Z_p <- as.vector(similarities %*% prob_cov) + 
    rep(sim_mean * (1 - sample_coverage), length(prob_cov))
  
  if (estimator %in% c("UnveilJ", "UnveilC", "UnveiliC")) {
    # concatenate the ordinariness of unobserved species, i.e.
    # their probability + sim_mean * the other probabilities 
    Z_p <- c(
      Z_p,
      prob_unv[-(1:s_obs)] + (1 - prob_unv[-(1:s_obs)]) * sim_mean
    )
    # Apply the naive estimator to the unveiled distribution
    if (q == 1) {
      # Limit value
      the_entropy <- as.numeric(-prob_unv %*% log(Z_p))
    } else {
      the_entropy <- as.numeric((1 - (prob_unv %*% (Z_p^(q - 1)))) / (q - 1))
    }
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
  
  ## MarconZhang ----
  if (estimator == "MarconZhang" | estimator == "Max") {
    V <- seq_len(sample_size - 1)
    # p_V_Ns is an array, containing (1 - (n_s-1)/(n-j)) for each species (lines) 
    # and all j from 1 to n-1
    p_V_Ns <- outer(abd, V, function(n_s, j) 1 - (n_s - 1) / (sample_size - j))
    # Useful values are products from j=1 to v, so prepare cumulative products
    p_V_Ns <- apply(p_V_Ns, 1, cumprod)
  }
  
  if (estimator == "ChaoShen" | estimator == "Max") {
    # Horvitz-Thomson multiplier
    HT_C_P <- prob_cov / (1 - (1 - prob_cov)^sample_size)
    # Force 0/0=0 and 0log0=0
    HT_C_P[prob_cov == 0] <- 0
    Z_p[Z_p == 0] <- 1
    ent_chao_shen <- (sum(HT_C_P * ln_q(1 / Z_p, q = q)))
  }
  
  if ((estimator == "MarconZhang" | estimator == "Max") & (q != 1)) {
    Zpqm1 <- Z_p^(q - 1)
    # Force 0^(q-1) = 0
    Zpqm1[Z_p == 0] <- 0
    K <- sum(prob_cov * Zpqm1)
    # Weights
    i <- seq_len(sample_size)
    w_vi <- (1 - sim_mean) * (i - q) / i
    w_v <- cumprod(w_vi)
    Taylor <- 1 + sum(
      prob * vapply(
        seq_len(s_obs), 
        FUN = S_v, 
        # Arguments
        abd = abd,
        sample_size = sample_size,
        w_v = w_v,
        p_V_Ns = p_V_Ns,
        # Value
        FUN.VALUE = 0
      )
    )
    FirstTerms <- prob_cov * (sim_mean + (1 - sim_mean) * prob_cov)^(q - 1)
    U <- Taylor - sum(FirstTerms)
    ent_marcon_zhang <- ((K + U - 1) / (1 - q))
  }
  
  if ((estimator == "MarconZhang" | estimator == "Max") & (q == 1)) {
    L <- -sum(prob_cov * log(Z_p))
    # Weights
    w_v <- ((1 - sim_mean)^V) / V
    Taylor <- 1 + sum(
      prob * vapply(
        seq_len(s_obs), 
        FUN = S_v, 
        # Arguments
        abd = abd,
        sample_size = sample_size,
        w_v = w_v,
        p_V_Ns = p_V_Ns,
        # Value
        FUN.VALUE = 0
      )
    )
    FirstTerms <- -prob_cov * log(sim_mean + (1 - sim_mean) * prob_cov)
    X <- Taylor - sum(FirstTerms)
    ent_marcon_zhang <- L + X
  }
  
  the_entropy <- switch (
    estimator,
    ChaoShen = ent_chao_shen,
    MarconZhang = ent_marcon_zhang,
    Max = max(ent_chao_shen, ent_marcon_zhang)
  )
  
  # Return
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


#' @rdname ent_similarity
#'
#' @export
ent_similarity.species_distribution <- function(
    x, 
    similarities = diag(sum(!(colnames(x) %in% non_species_columns))),
    q = 1, 
    estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang", 
                  "UnveilC", "UnveiliC", "naive"),
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Check species names
  is_species_column <- !(colnames(x) %in% non_species_columns)
  species_names <- colnames(x)[is_species_column]
  # Similarities may not be named
  if (is.null(colnames(similarities))) {
    if (ncol(similarities) != length(species_names)) {
      stop("If the similarity matrix is not named, then its size must fit the number of species.")
    }
  } else {
    if (length(setdiff(species_names, colnames(similarities))) != 0) {
      stop("Some species are missing in the similarity matrix.")    
    } else {
      # Filter and reorder the similarity matrix
      similarities <- similarities[species_names, species_names]
    }
  }
 
  if (gamma) {
    return(
      ent_gamma_similarity(
        species_distribution = x,
        similarities = similarities,
        q = q,
        estimator = estimator,
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        jack_alpha  = jack_alpha,
        jack_max = jack_max,
        as_numeric = FALSE
      )
    )
  } else {
    # Apply ent_similarity.numeric() to each site
    ent_similarity_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% non_species_columns)], 
      # Apply to each row
      MARGIN = 1,
      FUN = ent_similarity.numeric,
      # Arguments
      similarities = similarities,
      q = q,
      estimator = estimator,
      probability_estimator = probability_estimator,
      unveiling = unveiling,
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
        do.call(rbind.data.frame, ent_similarity_list)
      )
    )
  }
}


#' Sum of products weighted by w_v
#' 
#' Utility for Marcon-Zhang estimator
#'
#' @param species_index The species to consider (from 1 to `s_obs`)
#' @param abd An vector of abundances.
#' @param sample_size The sample size.
#' @param w_v A weight.
#' @param p_V_Ns An intermediate computation.
#'
#' @return A number.
#' @noRd
S_v <- function(
    species_index,
    abd,
    sample_size,
    w_v,
    p_V_Ns
    ) {
  v_used <- seq_len(sample_size - abd[species_index])
  return (sum(w_v[v_used] * p_V_Ns[v_used, species_index]))
}

#' Gamma entropy of a metacommunity
#' 
#' `species_distribution` is assumed to be a [species_distribution].
#' 
#' @return A tibble with the estimator used and the estimated entropy.
#' @noRd
ent_gamma_similarity <- function(
    species_distribution,
    similarities,
    q,
    estimator,
    probability_estimator,
    unveiling,
    jack_alpha,
    jack_max,
    as_numeric) {
  
  # Build the metacommunity
  abd <- metacommunity(
    species_distribution, 
    as_numeric = TRUE, 
    check_arguments = FALSE
  )
  if (is_integer_values(abd)) {
    # Sample coverage is useless
    sample_coverage <- NULL
  } else {
    # Non-integer values in the metacommunity. 
    # Calculate the sample coverage and change the estimator.
    sample_coverage <- coverage.numeric(
      colSums(
        species_distribution[
          , !(colnames(species_distribution) %in% non_species_columns)
        ]
      ),
      estimator = coverage_estimator,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
    if (!estimator %in% c("Marcon", "ChaoShen")) {
      estimator <- "Marcon"
    }
  }
  
  # Compute the entropy.
  the_entropy <- ent_similarity.numeric(
    abd,
    similarities = similarities,
    q = q,
    estimator = estimator,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    sample_coverage = sample_coverage,
    as_numeric = as_numeric,
    check_arguments = FALSE
  )
  # Add the site column
  if (!as_numeric) {
    the_entropy <- dplyr::bind_cols(
      site = "Metacommunity",
      the_entropy
    )
  }
  return(the_entropy)
}
