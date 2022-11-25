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
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
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
#' coverage(paracou_6_abd)
#' coverage_to_size(paracou_6_abd, sample_coverage = 0.9)
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
#' @param as_numeric If `TRUE`, a number is returned rather than a tibble.
#'
#' @export
coverage.numeric <- function(
    x, 
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    level = NULL, 
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  
  # Round values
  abd <- as.integer(round(x))
  # Eliminate zeros
  abd <- abd[abd > 0]
  # Calculate abundance distribution
  abd_distribution <- tapply(abd, abd, length)
  # singletons. Convert named number to number.
  s_1 <- as.numeric(abd_distribution["1"])
  sample_size <- sum(abd)
  
  # Coverage at the observed level ----
  if (is.null(level)) {
    # Estimate coverage at the observed level
    
    ##  No singletons: coverage=1 ----
    if (is.na(s_1)) {
      if (as_numeric) {
        return(1)
      } else {
        return(
          tibble::tibble_row(
            estimator = "No singleton", 
            coverage = 1
          )
        )  
      }
    }
    
    ## Singletons only: coverage=0 ----
    if (s_1 == sample_size) {
      warning ("Sample coverage is 0, most bias corrections will return NaN.")
      if (as_numeric) {
        return(0)
      } else {
        return(
          tibble::tibble_row(
            estimator = "Singletons only", 
            coverage = 0
          )
        )  
      }
    }
  
    ## Zhang & Huang's estimator ----
    if (estimator == "ZhangHuang") {
      prob <- abd/sample_size
      if (any(prob >= .5)) {
        warning ("Zhang-Huang's sample coverage cannot be estimated because one probability is over 1/2. Chao's estimator is returned.")
        estimator <- "Chao"
      } else {
        nu <- as.integer(names(abd_distribution))
        # Use nu %% 2 * 2 - 1 for (-1)^(Nu + 1)
        sample_coverage <- 1 - sum((nu %% 2 * 2 - 1) / choose(sample_size, nu) * abd_distribution)
        if (as_numeric) {
          return(sample_coverage)
        } else {
          return(
            tibble::tibble_row(
              estimator = estimator, 
              coverage = sample_coverage
            )
          )  
        }
      }    
    }
    
    ## Chao's estimator ----
    if (estimator == "Chao") {
      sample_coverage <- 1 - s_1 / sample_size * (1-chao_A(abd))
      if (as_numeric) {
        return(sample_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            coverage = sample_coverage
          )
        )  
      }
    }
    
    ## Turing's estimator ----
    if (estimator == "Turing") {
      sample_coverage <- 1 - s_1 / sample_size
      if (as_numeric) {
        return(sample_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            coverage = sample_coverage
          )
        )  
      }
    }
  }
  
  # Coverage at a chosen level ----
  if (!is.null(level)) {
    # level must be an integer. check_divent_args() may have accepted a value between 0 and 1
    if (level <=1) stop("level must be an integer > 1.")

    ## Good's estimator ----
    if (estimator == "Good") {
      if (level >= sample_size) stop("Good's estimator only allows interpolation: level must be less than the observed community size.")
      sample_coverage <- 1 - EntropyEstimation::GenSimp.z(abd, level)
      if (as_numeric) {
        return(sample_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            coverage = sample_coverage
          )
        )  
      }
    }
    
    ## Chao's estimator ----
    if (estimator == "Chao") {
      if (level < sample_size) {
        ### Interpolation ----
        abd_restricted <- abd[(sample_size - abd) >= level]
        sample_coverage <- 1 - sum(
          abd_restricted / sample_size * 
          exp(lgamma(sample_size - abd_restricted + 1) - 
          lgamma(sample_size - abd_restricted - level + 1) - 
          lgamma(sample_size) + lgamma(sample_size - level))
        )
      } else {
        ### Extrapolation ----
        if (is.na(s_1)) {
          # No singletons, coverage=1
          if (as_numeric) {
            return(1)
          } else {
            return(
              tibble::tibble_row(
                estimator = "No singleton", 
                coverage = 1
              )
            )  
          }
        } else {
          sample_coverage <- 1 - 
            s_1 / sample_size * 
            (1 - chao_A(abd))^(level - sample_size + 1)
        }
      }
      if (as_numeric) {
        return(sample_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            coverage = sample_coverage
          )
        )  
      }
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
#' @export
coverage_to_size <- function(x, ...) {
  UseMethod("coverage_to_size")
}

#' @rdname sample_coverage
#'
#' @param sample_coverage The target sample coverage.
#'
#' @export
coverage_to_size.numeric <- function(
    x, 
    sample_coverage,
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  
  # Round values
  abd <- as.integer(round(x))
  # Eliminate zeros
  abd <- abd[abd > 0]
  # Calculate abundance distribution
  abd_distribution <- tapply(abd, abd, length)
  # Singletons. Convert named number to number.
  if (is.na(abd_distribution["1"])) {
    s_1 <- 0 
  } else {
    s_1 <- as.numeric(abd_distribution["1"])
  }
  sample_size <- sum(abd)
  
  # Singletons only
  if (s_1 == sample_size) {
    stop("Sample coverage is 0.")
  }
  
  # Actual coverage
  sample_coverage_actual <- coverage.numeric(abd, check_arguments = FALSE)$coverage
  
  if (sample_coverage >= sample_coverage_actual) {
    # Extrapolation
    size <- round(
      sample_size + 
        (log(sample_size / s_1) + log(1 - sample_coverage)) / 
          log(1 - chao_A(abd)
      ) - 1
    )
  } else {
    # Interpolation. Numeric resolution: minimize the function delta
    size <- round(
      stats::optimize(
        chao_delta,
        abd = abd,
        target_coverage = sample_coverage, 
        lower = 1, 
        upper = sample_size
      )$minimum
    )
  }
  return(tibble::tibble_row(sample_coverage = sample_coverage, size = size))
  
}


#' @rdname sample_coverage
#' 
#' @export
coverage_to_size.abundances <- function(
    x, 
    sample_coverage,
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()

  # Apply coverage_to_size.numeric() to each site
  size_list <- apply(
    # Eliminate site and weight columns
    x[, !(colnames(x) %in% c("site", "weight"))], 
    # Apply to each row
    MARGIN = 1,
    FUN = coverage_to_size.numeric,
    # Arguments
    sample_coverage = sample_coverage,
    check_arguments = FALSE
  )
  
  return(
    # Make a tibble with site, estimator and sample-coverage
    tibble::tibble(
      # Do not assume column site exists
      x[colnames(x) == "site"],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, size_list)
    )
  )
}


#' Chao's A
#' 
#' Helper for Chao's estimators.
#' 
#' A's formula \insertCite{@Chao2015@, eq. 6b}{divent}) depends on the presence of singletons and doubletons.
#' 
#' @noRd
#'
#' @param abd A vector of positive integers (not checked).
#'
#' @return The value of A.
chao_A <- function(abd) {
  
  # Calculate abundance distribution
  abd_distribution <- tapply(abd, abd, length)
  s_1 <- as.numeric(abd_distribution["1"])
  s_2 <- as.numeric(abd_distribution["2"])
  sample_size <- sum(abd)
  
  # Calculate A
  if (is.na(s_1)) {
    A <- 0
  } else {
    # Use Chao1 estimator to evaluate the number of unobserved species
    if (is.na(s_2)) {
      s_0 <- (sample_size - 1) * s_1 * (s_1 - 1) / 2 / sample_size
    } else {
      s_0 <- (sample_size - 1) * s_1^2 / 2 / sample_size / s_2
    }
    A <- 1- sample_size * s_0 / (sample_size * s_0 + s_1)
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
chao_delta <- function(abd, size, target_coverage) {
  abs(
    coverage.numeric(
      abd, 
      estimator = "Chao", 
      level = size, 
      check_arguments = FALSE
    )$coverage - target_coverage
  )
}