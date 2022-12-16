#' Probabilities of Species
#' 
#' Estimate actual probabilities of species from a sample
#' 
#' `probabilities()` estimates a probability distribution from a sample.
#' If the `estimator` is not "naive", the observed abundance distribution is used 
#' to estimate the actual probability distribution. 
#' The list of species is changed:
#' zero-abundance species are cleared, and some unobserved species are added. 
#' First, observed species probabilities are estimated following 
#' \insertCite{Chao2003;textual}{divent}, i.e. input probabilities are multiplied by 
#' the sample coverage, or according to more sophisticated models: 
#' \insertCite{Chao2013;textual}{divent}, a single-parameter model, or 
#' \insertCite{Chao2015;textual}{divent}, a two-parameter model. 
#' The total probability of observed species equals the sample coverage. 
#' Then, the distribution of unobserved species can be unveiled: their number 
#' is estimated according to the `richness_estimator`. 
#' The coverage deficit (1 minus the sample coverage) is shared by the unobserved 
#' species equally (`unveiling` = "uniform", \insertCite{Chao2013}{divent}) or 
#' according to a geometric distribution (`unveiling` = "geometric", \insertCite{Chao2015}{divent}).
#'  
#'
#' @inheritParams check_divent_args
#' @param x An object. It may be:
#' 
#' - a numeric vector containing abundances. It may be named to track species names.
#' - an object of class [species_distribution].
#' 
#' @param ... Unused.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Just transform abundances into probabilities
#' probabilities(paracou_6_abd)
#' # Estimate the distribution of probabilities from observed abundances (unveiled probabilities)
#' prob_unv <- probabilities(
#'   paracou_6_abd, 
#'   estimator = "Chao2015", 
#'   unveiling = "geometric",
#'   richness_estimator = "jackknife"
#' )
#' # Estimate entropy from the unveiled probabilities
#' ent_shannon(prob_unv)
#' # Identical to
#' ent_shannon(paracou_6_abd, estimator = "UnveilJ")
#' 
#' @name probabilities
NULL


#  Probabilities ----

#' @rdname probabilities
#'
#' @export
probabilities <- function(x, ...) {
  UseMethod("probabilities")
}


#' @rdname probabilities
#'
#' @param estimator One of the estimators of a probability distribution: 
#' "naive" (the default value), or "Chao2013", "Chao2015", "ChaoShen" to estimate
#' the probabilities of the observed species in the asymptotic distribution.
#' @param unveiling One of the possible unveiling methods to estimate the probabilities 
#' of the unobserved species: "none" (default, no species is added), "uniform" 
#' (all unobserved species have the same probability) or "geometric" (the 
#' unobserved species distribution is geometric).
#' @param richness_estimator An estimator of richness to evaluate the total number 
#' of species. "jackknife" is the default value. 
#' An alternative is "rarefy" to estimate the number of species such that the 
#' entropy of the asymptotic distribution rarefied to the observed sample size equals
#' the actual entropy of the data.
#' @param coverage_estimator An estimator of sample coverage used by [coverage]:
#' "ZhangHuang" (the default value), \code{"Chao"}, \code{"Turing"} or \code{"Good"}.
#' @param q The order of diversity. Default is 0 for richness. 
#' Used only to estimate asymptotic probability distributions when argument 
#' `richness_estimator` is "rarefy". Then, the number of unobserved species is 
#' fitted so that the entropy of order q of the asymptotic probability distribution 
#' at the observed sample size equals the actual entropy of the data.
#'
#' @export
probabilities.numeric <- function(
    x, 
    estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"), 
    jack_alpha = 0.05, 
    jack_max = 10, 
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q = 0,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {
  
  # Check the data ----
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  
  
  # Save or generate the species names
  species_names <- names(x)

  # Naive estimator ----
  if (estimator == "naive") {
    # Just normalize so that x sums to 1
    the_prob <- x / sum(x)
  # Other estimators ----
  } else {
    # Integer abundances are required by all non-naive estimators
    if (!is_integer_values(x)) {
      warning(
        "Integer abundance values are required to estimate community probabilities. Abundances have been rounded."
      )
    }
    abd_int <- round(x)
    
    # Eliminate 0 and calculate elementary statistics
    abd <- abd_int[abd_int > 0]
    s_obs <- length(abd)
    sample_size <- sum(abd)
    prob <- abd / sample_size
    # Sample coverage
    sample_coverage <- coverage.numeric(
      abd, 
      estimator = coverage_estimator,
      as_numeric = TRUE,
      check_arguments = FALSE
      )
    if (
      estimator == "Chao2015" | 
      unveiling != "none" | 
      richness_estimator == "Rarefy") {
      # Sample coverage of order 2 required
      s_1 <- sum(abd == 1)
      s_2 <- sum(abd == 2)
      if (s_2 == 0) {
        s_1 <- max(s_1 - 1, 0)
        s_2 <- 1
      }
      s_3 <- max(sum(abd == 3), 1)
      # 1 minus sample coverage (i.e. Coverage Deficit) of order 2
      coverage_deficit_2 <-
        s_2 / choose(sample_size, 2) * 
        ((sample_size - 2) * s_2 / ((sample_size - 2) * s_2 + 3 * s_3))^2
    }
    
    ## Tune the probabilities of observed species ----
    if (sample_coverage == 0 | sample_coverage == 1) {
      # Sample coverage equal to 1, do not tune. If 0, unable to tune.
      prob_tuned <- prob
    } else {
      prob_tuned <- NA
      if (estimator == "ChaoShen") {
        prob_tuned <- sample_coverage * prob
      }
      if (estimator == "Chao2013") {
        # Single parameter estimation, Chao et al. (2013)
        denominator <- sum(prob * (1 - prob)^sample_size)
        if (denominator == 0) {
          # sample_size is too big so denominator equals 0. 
          # Just multiply by coverage.
          prob_tuned <- sample_coverage * prob
        } else {
          # General case
          lambda <- (1 - sample_coverage) / denominator
          prob_tuned <- prob * (1 - lambda * (1 -prob)^sample_size)
        }      
      } 
      if (estimator == "Chao2015")  {
        # Two parameters, Chao et al. (2015). 
        # Estimate theta. Set it to 1 if impossible.
        theta <- tryCatch(
          stats::optimize(
            theta_solve, 
            interval = c(0, 1), 
            prob = prob, 
            abd = abd, 
            sample_size = sample_size, 
            sample_coverage = sample_coverage, 
            coverage_deficit_2 = coverage_deficit_2
          )$min, 
          error = function(e) {1}
        )
        lambda <- (1 - sample_coverage) / sum(prob * exp(-theta * abd))
        prob_tuned <- prob * (1 - lambda * exp(-theta * abd))
      }
    }
    # Restore the species names
    if (is.null(species_names)) {
      # No names: create them
      names(prob_tuned) <- paste("sp", seq_along(prob_tuned), sep = "")
    } else {
      # Restore the names
      names(prob_tuned) <- species_names[abd_int > 0]
    }
    
    
    ## Estimate the number of unobserved species ----
    if (richness_estimator == "rarefy") {
      if (unveiling == "none") {
        stop("Arguments richness_estimator='rarefy' and unveiling='none' are not compatible")
      }
      # Estimation of the number of unobserved species to initialize optimization
      s_0 <- div_richness.numeric(
        abd, 
        estimator = "jackknife", 
        as_numeric = TRUE, 
        check_arguments = FALSE
      ) - s_obs
      # Estimate the number of unobserved species by iterations
      ent_target <- ent_tsallis.numeric(
        abd, 
        q = q, 
        estimator = "naive", 
        as_numeric = TRUE, 
        check_arguments = FALSE
      )
      s_0 <- round(
        tryCatch(
          stats::optimize(
            rarefaction_bias, 
            interval = c(0, 2 * s_0), 
            abd = abd, 
            prob_tuned = prob_tuned, 
            sample_coverage = sample_coverage, 
            coverage_deficit_2 = coverage_deficit_2, 
            q = q, 
            unveiling = unveiling, 
            ent_target = ent_target
          )$minimum,
          error = function(e) {s_0}
        )
      )
    } else {
      s_est <- ceiling(
        div_richness.numeric(
          abd, 
          estimator = richness_estimator, 
          jack_alpha = jack_alpha, 
          jack_max = jack_max,
          as_numeric = TRUE, 
          check_arguments = FALSE
        )
      )
      s_0 <- s_est - s_obs
    }
    
    ## Distribution of unobserved species ----
    if (s_0) {
      if (unveiling == "none") {
        the_prob <- prob_tuned
      } else {
        the_prob <- c(
          prob_tuned, 
          estimate_prob_s_0(
            unveiling = unveiling, 
            prob_tuned = prob_tuned, 
            s_0 = s_0, 
            sample_coverage = sample_coverage, 
            coverage_deficit_2 = coverage_deficit_2
          )
        )
      }
    } else {
      the_prob <- prob_tuned
    }
  }
  
  # Set the class and return ----
  if (as_numeric) {
    return(the_prob)
  } else {
    the_probabilities <- as_species_distribution(the_prob)
    class(the_probabilities) <- c("probabilities", class(the_probabilities))
    return(the_probabilities)
  }
}


#' @rdname probabilities
#'
#' @export
probabilities.abundances <- function(
    x, 
    estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"), 
    jack_alpha = 0.05, 
    jack_max = 10, 
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q = 0,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")

  # Apply probabilities.numeric() to each site
  probabilities_list <- apply(
    # Eliminate site and weight columns
    x[, !(colnames(x) %in% c("site", "weight"))], 
    # Apply to each row
    MARGIN = 1,
    FUN = probabilities.numeric,
    # Arguments
    estimator = estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha = jack_alpha, 
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    q = q,
    check_arguments = FALSE
  )
  
  # Bind the rows
  the_probabilities <- dplyr::bind_rows(probabilities_list)
  # Restore the site names
  if (("site" %in% colnames(x))) the_probabilities$site <- x$site
  # Replace NA's due to binding by zeros
  the_probabilities <- dplyr::mutate(
    the_probabilities,
    dplyr::across(
      .cols = dplyr::everything(),
      .fns = ~ ifelse(is.na(.x), 0, .x))
  )
  return(the_probabilities)
}


#  Utilities ----

#' Solve the theta parameter of Chao et al. (2015)
#' 
#' Utilities for [probabilities.numeric].
#' 
#' Code inspired from JADE function DetAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
#' 
#' @noRd
#'
#' @param theta The parameter to solve.
#' @param prob A vector of probabilities (not checked).
#' @param abd A vector of positive integers (not checked).
#' @param sample_size The number of individuals in the sample.
#' @param sample_coverage The sample coverage. 
#' @param coverage_deficit_2 The coverage deficit of order 2.
#'
#' @return The value of the parameter theta to minimize.
theta_solve <- function(
    theta, 
    prob, 
    abd, 
    sample_size, 
    sample_coverage, 
    coverage_deficit_2) {
  
  lambda <- (1 - sample_coverage) / sum(prob * exp(-theta * abd))
  return(
    abs(
      sum((prob * (1 - lambda * exp(-theta * abd)))^2) - 
        sum(choose(abd, 2) / choose(sample_size, 2)) + coverage_deficit_2
    )
  )
}


#' Solve the beta parameter of Chao et al. (2015)
#'
#' Utilities for [probabilities.numeric].
#' 
#' Code inspired from JADE function UndAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
#' 
#' @noRd
#' 
#' @param beta The parameter to solve.
#' @param r The squared coverage deficit divided by the coverage deficit of order 2.
#' @param i The sequence from 1 to the number of species.
#'
#' @return The value of the parameter beta to minimize.
beta_solve <- function(beta, r, i) {
  return(abs(sum(beta^i)^2 / sum((beta^i)^2) - r))
}


#' Unobserved Species Distribution
#'
#' Utilities for [probabilities.numeric].
#' 
#' @noRd
#' 
#' @param unveiling The unveiling method.
#' @param prob_tuned The tuned distribution of probabilities.
#' @param s_0 The number of unobserved species.
#' @param sample_coverage The sample coverage. 
#' @param coverage_deficit_2 The coverage deficit of order 2.
#'
#' @return The distribution of probabilities of unobserved species.
estimate_prob_s_0 <- function(
    unveiling, 
    prob_tuned, 
    s_0, 
    sample_coverage, 
    coverage_deficit_2) {

  the_prob_s_0 <- NA
  if (unveiling == "geometric") {
    if (s_0 == 1) {
      # A single unobserved species
      the_prob_s_0 <- 1 - sample_coverage
    } else {
      r <- (1 - sample_coverage)^2 / coverage_deficit_2
      i <- seq_len(s_0)
      beta <-  tryCatch(
        stats::optimize(
          beta_solve, 
          lower = (r - 1) / (r + 1), 
          upper = 1, 
          tol = .Machine$double.eps, 
          r, 
          i
        )$min, 
        error = function(e) {(r - 1) / (r + 1)}
      )
      alpha <- (1 - sample_coverage) / sum(beta^i)
      the_prob_s_0 <- alpha * beta^i
      # Sometimes fails when the distribution is very uneven (sometimes r < 1) 
      # Then, fall back to the uniform distribution
      if (any(is.na(the_prob_s_0)) | any(the_prob_s_0 <= 0)) {
        unveiling <- "uniform"
      }
    }
  }      
  if (unveiling == "uniform") {
    # Add s_0 unobserved species with equal probabilities
    the_prob_s_0 <- rep((1 - sum(prob_tuned)) / s_0, s_0)
  }
  if (any(is.na(the_prob_s_0))) {
    warning("Unveiling method was not recognized")
    return(NA)
  } else {
    names(the_prob_s_0) <- paste("Unobs_sp", seq_along(the_prob_s_0), sep = "_")
    return(the_prob_s_0)
  }         
}


#' Rarefaction Bias
#' 
#' Departure of the rarefied entropy from the target entropy
#'
#' Utilities for [probabilities.numeric].
#'
#' @noRd
#'
#' @param s_0 The number of unobserved species.
#' @param abd The abundances of species.
#' @param prob_tuned The tuned distribution of probabilities.
#' @param sample_coverage The sample coverage.
#' @param coverage_deficit_2 The coverage deficit of order 2.
#' @param q The order of entropy to fit.
#' @param unveiling The unveiling method.
#' @param ent_target Target entropy.
#'
#' @return The departure of the rarefied entropy from the target entropy.
rarefaction_bias <- function(
    s_0,
    abd, 
    prob_tuned, 
    sample_coverage, 
    coverage_deficit_2, 
    q, 
    unveiling, 
    ent_target) {
  
  abd <- abd[abd > 0]
  sample_size <- sum(abd)
  # Unobserved species
  prob_s_0 <- estimate_prob_s_0(
    unveiling, 
    prob_tuned, 
    s_0, 
    sample_coverage, 
    coverage_deficit_2
  )
  # Full distribution of probabilities
  prob <- c(prob_tuned, prob_s_0)
  # abundances_freq_count at level = sample_size
  s_nu <- vapply(
    seq_len(sample_size), 
    function(nu) {
      sum(
        exp(
          lchoose(sample_size, nu) + nu * log(prob) + 
            (sample_size - nu) * log(1 - prob)
        )
      )
    }, 
    FUN.VALUE=0.0
  )
  # Get entropy at level=sample_size and calculate the bias
  if (q == 1) {
    the_ent_bias <- abs(
      sum(
        -seq_len(sample_size) / sample_size * 
          log(seq_len(sample_size) / sample_size) * s_nu
      ) 
      - ent_target
    )
  } else {
    the_ent_bias <- abs(
      (sum((seq_len(sample_size)/sample_size)^q * s_nu) - 1) / (1 - q) 
      - ent_target
    )
  }
  return(the_ent_bias)
}
