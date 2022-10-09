#' Sample Coverage of a Community
#' 
#' `coverage()` calculates an estimator of the sample coverage of a community 
#' described by its abundance vector. 
#' `coverage_to_size()` estimates the sample size corresponding to the chosen 
#' sample coverage.
#'
#' @param x An object.
#' @param abundances An object of class `abundances`.
#' @param ... Unused.
#' @param check_arguments If `TRUE`, check the arguments of the function. 
#' 
#' @return `coverage()` returns a named number equal to the calculated sample coverage.
#' The name is that of the estimator used. 
#' `coverage_to_size()` returns a number equal to the sample size corresponding
#' to the chosen sample coverage.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' coverage(asm_paracou_6)
#' 
#' @name sample_coverage
NULL


#' @rdname sample_coverage
#'
#' @export
coverage <- function(x, ...) {
  UseMethod("coverage")
}

#' @rdname sample_coverage
#' 
#' @param estimator A string containing one of the possible estimators: 
#' "ZhangHuang", "Chao", "Turing", "Good".
#' @param level The level of interpolation or extrapolation, i.e. an abundance.
#'
#' @export
coverage.numeric <- function(
    x, 
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    level = NULL, 
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  
  # Round values
  abundances <- as.integer(round(x))
  # Eliminate zeros
  abundances <- abundances[abundances > 0]
  # Calculate abundance distribution
  abundance_distribution <- tapply(abundances, abundances, length)
  # singletons. Convert named number to number.
  singletons <- as.numeric(abundance_distribution["1"])
  sample_size <- sum(abundances)
  
  if (is.null(level)) {
    # Estimate coverage at the observed level
    
    # No singletons: coverage=1
    if (is.na(singletons)) {
      return(tibble::tibble_row(estimator = "No singleton", coverage = 1))
    }
    
    # singletons only
    if (singletons == sample_size) {
      warning ("Sample coverage is 0, most bias corrections will return NaN.")
      return(tibble::tibble_row(estimator = "Singletons only", coverage = 0))
    }
  
    if (estimator == "ZhangHuang") {
      probabilities <- abundances/sample_size
      if (any(probabilities >= .5)) {
        warning ("Zhang-Huang's sample coverage cannot be estimated because one probability is over 1/2. Chao's estimator is returned.")
        estimator <- "Chao"
      } else {
        nu <- as.integer(names(abundance_distribution))
        # Use nu %% 2 * 2 - 1 for (-1)^(Nu + 1)
        sample_coverage <- 1 - sum((nu %% 2 * 2 - 1) / choose(sample_size, nu) * abundance_distribution)
        return(tibble::tibble_row(estimator, coverage = sample_coverage))
      }    
    }
    if (estimator == "Chao") {
      sample_coverage <- 1 - singletons / sample_size * (1-chao_A(abundances))
      return(tibble::tibble_row(estimator, coverage = sample_coverage))
    }
    if (estimator == "Turing") {
      sample_coverage <- 1 - singletons / sample_size
      return(tibble::tibble_row(estimator, coverage = sample_coverage))
    }
    
  } else {
    # Chose level. Must be an integer. check_divent_args() may have accepted a value between 0 and 1
    if (level <=1) stop("level must be an integer > 1.")

    if (estimator == "Best") estimator <- "Chao"
    # Extrapolation allowed.
    
    if (estimator == "Good") {
      if (level >= sample_size) stop("The Good estimator only allows interpolation: level must be less than the observed community size.")
      sample_coverage <- 1 - EntropyEstimation::Geabundancesimp.z(abundances, level)
      return(tibble::tibble_row(estimator, coverage = sample_coverage))
    }
    if (estimator == "Chao") {
      if (level < sample_size) {
        # Interpolation
        abundancesRestricted <- abundances[(sample_size - abundances) >= level]
        sample_coverage <- 1 - sum(
          abundancesRestricted/sample_size * exp(lgamma(sample_size - abundancesRestricted + 1) - 
          lgamma(sample_size - abundancesRestricted - level + 1) - 
          lgamma(sample_size) + lgamma(sample_size - level))
        )
      } else {
        # Extrapolation
        if (is.na(singletons)) {
          # No singletons, C=1
          return(tibble::tibble_row(estimator = "No singleton", coverage = 1))
        } else {
          sample_coverage <- 1 - singletons / sample_size * (1 - chao_A(abundances))^(level - sample_size + 1)
        }
      }
      return(tibble::tibble_row(estimator, coverage = sample_coverage))
    }
  }
}


#' @rdname sample_coverage
#' 
#' @export
coverage.abundances <- function(
    x, 
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    level = NULL, 
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  
  # Apply coverage.numeric() to each site
  coverage_list <- apply(
    # Eliminate site and weight columns
    x[, !(colnames(x) %in% c("site", "weight"))], 
    # Apply to each row
    MARGIN = 1,
    FUN = coverage.numeric,
    # Arguments
    estimator = estimator,
    level = level,
    check_arguments = FALSE
  )
  
  return(
    # Make a tibble with site, estimator and sample-coverage
    tibble::tibble(
      # Do not assume column site exists
      x[colnames(x) == "site"],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, coverage_list)
    )
  )
}
  
  
#' @rdname sample_coverage
#'
#' @param sample_coverage The target sample coverage.
#'
#' @export
coverage_to_size <- function(
    abundances, 
    sample_coverage, 
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  
  # Round values
  abundances <- as.integer(round(abundances))
  # Eliminate zeros
  abundances <- abundances[abundances > 0]
  # Calculate abundance distribution
  abundance_distribution <- tapply(abundances, abundances, length)
  # singletons. Convert named number to number.
  if (is.na(abundance_distribution["1"])) {
    singletons <- 0 
  } else {
    singletons <- as.numeric(abundance_distribution["1"])
  }
  sample_size <- sum(abundances)
  
  # singletons only
  if (singletons == sample_size) {
    stop("Sample coverage is 0.")
  }
  
  # Actual coverage
  sample_coverage_actual <- coverage(abundances, check_arguments = FALSE)
  
  if (sample_coverage >= sample_coverage_actual) {
    # Extrapolation
    size <- round(
      sample_size + 
        (log(sample_size / singletons) + log(1 - sample_coverage)) / 
          log(1 - chao_A(abundances)
      ) - 1
    )
  } else {
    # Interpolation. Numeric resolution: minimize the function delta
    size <- round(
      stats::optimize(
        chao_delta, 
        target_coverage = sample_coverage, 
        lower = 1, 
        upper = sample_size
      )$minimum
    )
  }
  return(size)
}


#' Chao's A
#' 
#' Helper for Chao's estimators.
#' 
#' A's formula depends on the presence of singletons and doubletons.
#' 
#' @noRd
#'
#' @param abundances A vector of positive integers (not checked).
#'
#' @return The value of A.
chao_A <- function(abundances) {
  
  # Calculate abundance distribution
  abundance_distribution <- tapply(abundances, abundances, length)
  singletons <- as.numeric(abundance_distribution["1"])
  doubletons <- as.numeric(abundance_distribution["2"])
  sample_size <- sum(abundances)
  
  # Calculate A
  if (is.na(singletons)) {
    A <- 0
  } else {
    # Use Chao1 estimator to evaluate the number of unobserved species
    if (is.na(doubletons)) {
      species_unobserved <- (sample_size - 1) * singletons * (singletons - 1) / 2 / sample_size
    } else {
      species_unobserved <- (sample_size - 1) * singletons^2 / 2 / sample_size / doubletons
    }
    A <- 1- sample_size * species_unobserved / (sample_size * species_unobserved + singletons)
  }
  
  return(A)
}


#' Departure of actual sample coverage from target coverage
#' 
#' Helper for `coverage_2_size()`
#' 
#' @noRd
#'
#' @param size The size of the sample. Adjusted to minimize `delta()`.
#' @param target_coverage The sample coverage to reach by adjusting size.
#'
#' @return The departure of actual sample coverage from target coverage.
chao_delta <- function(size, target_coverage) {
  abs(
    coverage(
      abundances, 
      estimator = "Chao", 
      level = size, 
      check_arguments = FALSE
    ) - target_coverage
  )
}
