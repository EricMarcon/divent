#' Sample Coverage of a Community
#'
#' `coverage()` calculates an estimator of the sample coverage of a community
#' described by its abundance vector.
#' `coverage_to_size()` estimates the sample size corresponding to the chosen
#' sample coverage.
#'
#' The sample coverage \eqn{C} of a community is the total probability of
#' occurrence of the species observed in the sample.
#' \eqn{1-C} is the probability for an individual of the whole community to
#' belong to a species that has not been sampled.
#'
#' The historical estimator is due to Turing \insertCite{Good1953}{divent}.
#' It only relies on singletons (species observed only once).
#' Chao's \insertCite{Chao2010a}{divent} estimator uses doubletons too and
#' Zhang-Huang's \insertCite{Chao1988,Zhang2007}{divent} uses the whole
#' distribution.
#'
#' If `level` is not `NULL`, the sample coverage is interpolated or extrapolated.
#' Interpolation by the Good estimator relies on the equality between sampling
#' deficit and the generalized Simpson entropy \insertCite{Good1953}{divent}.
#' The \insertCite{Chao2014;textual}{divent} estimator allows extrapolation,
#' reliable up a level equal to the double size of the sample.
#'
#' @inheritParams check_divent_args
#' @param x An object.
#' @param ... Unused.
#'
#' @returns `coverage()` returns a named number equal to the calculated sample coverage.
#' The name is that of the estimator used.
#'
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
#' @name coverage
NULL


#' @rdname coverage
#'
#' @export
coverage <- function(x, ...) {
  UseMethod("coverage")
}

#' @rdname coverage
#'
#' @param estimator An estimator of the sample coverage.
#' "ZhangHuang" is the most accurate but does not allow choosing a `level`.
#' "Good"'s estimator only allows interpolation, i.e. estimation of the coverage
#' of a subsample.
#' "Chao" allows estimation at any `level`, including extrapolation.
#' "Turing" is the simplest estimator, equal to 1 minus the number of singletons
#' divided by the sample size.
#' @param level The level of interpolation or extrapolation, i.e. an abundance.
#'
#' @export
coverage.numeric <- function(
    x,
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    level = NULL,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) check_divent_args()

  # Round values
  abd <- as.integer(round(x))
  # Eliminate zeros
  abd <- abd[abd > 0]
  # Calculate abundance distribution
  abd_distribution <- tapply(abd, INDEX = abd, FUN = length)
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
      cli::cli_alert_warning(
        "Sample coverage is 0, most estimators will return {.code NaN}."
      )
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
        estimator <- "Chao"
      } else {
        nu <- as.integer(names(abd_distribution))
        # Use nu %% 2 * 2 - 1 for (-1)^(Nu + 1)
        the_coverage <- 1 - sum((nu %% 2 * 2 - 1) / choose(sample_size, nu) * abd_distribution)
        if (as_numeric) {
          return(the_coverage)
        } else {
          return(
            tibble::tibble_row(
              estimator = estimator,
              coverage = the_coverage
            )
          )
        }
      }
    }

    ## Chao's estimator ----
    if (estimator == "Chao") {
      the_coverage <- 1 - s_1 / sample_size * (1 - chao_A(abd))
      if (as_numeric) {
        return(the_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator,
            coverage = the_coverage
          )
        )
      }
    }

    ## Turing's estimator ----
    if (estimator == "Turing" | estimator == "Good") {
      the_coverage <- 1 - s_1 / sample_size
      if (as_numeric) {
        return(the_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator,
            coverage = the_coverage
          )
        )
      }
    }
  }

  # Coverage at a chosen level ----
  if (!is.null(level)) {
    # level must be an integer. check_divent_args() may have accepted a value between 0 and 1
    if (level <= 1) stop("level must be an integer > 1.")

    ## Turing or Zhang & Huang's estimator ----
    if (estimator == "Turing" | estimator == "ZhangHuang") {
      cli::cli_alert_warning(
        "Turing and ZhangHuang estimators do not allow interpolation or extrapolation."
      )
      cli::cli_alert("Chao's estimator is used.")
      estimator <- "Chao"
    }

    ## Good's estimator ----
    if (estimator == "Good") {
      if (level >= sample_size) {
        cli::cli_alert_warning("Good's estimator only allows interpolation.")
        cli::cli_alert("Chao's estimator is used.")
        estimator <- "Chao"
      } else {
        the_coverage <- 1 - EntropyEstimation::GenSimp.z(abd, level)
        if (as_numeric) {
          return(the_coverage)
        } else {
          return(
            tibble::tibble_row(
              estimator = estimator,
              level = level,
              coverage = the_coverage
            )
          )
        }
      }
    }

    ## Chao's estimator ----
    if (estimator == "Chao") {
      if (level < sample_size) {
        ### Interpolation ----
        abd_restricted <- abd[(sample_size - abd) >= level]
        the_coverage <- 1 - sum(
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
                level = level,
                coverage = 1
              )
            )
          }
        } else {
          the_coverage <- 1 -
            s_1 / sample_size *
            (1 - chao_A(abd))^(level - sample_size + 1)
        }
      }
      if (as_numeric) {
        return(the_coverage)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator,
            level = level,
            coverage = the_coverage
          )
        )
      }
    }
  }
}


#' @rdname coverage
#'
#' @export
coverage.abundances <- function(
    x,
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    level = NULL,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) check_divent_args()

  # Apply coverage.numeric() to each site
  coverage_list <- apply(
    # Eliminate site and weight columns
    x[, !colnames(x) %in% non_species_columns],
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
      # Restore non-species columns
      x[colnames(x) %in% non_species_columns],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, coverage_list)
    )
  )
}


#' @rdname coverage
#'
#' @export
coverage_to_size <- function(x, ...) {
  UseMethod("coverage_to_size")
}

#' @rdname coverage
#'
#' @param sample_coverage The target sample coverage.
#'
#' @export
coverage_to_size.numeric <- function(
    x,
    sample_coverage,
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) check_divent_args()

  # Round values
  abd <- as.integer(round(x))
  # Eliminate zeros
  abd <- abd[abd > 0]
  # Calculate abundance distribution
  abd_distribution <- tapply(abd, INDEX = abd, FUN = length)
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
  sample_coverage_actual <- coverage.numeric(
    abd,
    estimator = estimator,
    as_numeric = TRUE,
    check_arguments = FALSE
  )

  if (sample_coverage >= sample_coverage_actual) {
    # Extrapolation
    the_size <- round(
      sample_size +
        (log(sample_size / s_1) + log(1 - sample_coverage)) /
          log(1 - chao_A(abd)
      ) - 1
    )
  } else {
    # Interpolation. Numeric resolution: minimize the function delta
    the_size <- round(
      stats::optimize(
        chao_delta,
        abd = abd,
        target_coverage = sample_coverage,
        lower = 1,
        upper = sample_size
      )$minimum
    )
  }

  if (as_numeric) {
    return(the_size)
  } else {
    return(
      tibble::tibble_row(
        sample_coverage = sample_coverage,
        size = the_size
      )
    )
  }
}


#' @rdname coverage
#'
#' @export
coverage_to_size.abundances <- function(
    x,
    sample_coverage,
    estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) check_divent_args()

  # Apply coverage_to_size.numeric() to each site
  size_list <- apply(
    # Eliminate site and weight columns
    x[, !colnames(x) %in% non_species_columns],
    # Apply to each row
    MARGIN = 1,
    FUN = coverage_to_size.numeric,
    # Arguments
    sample_coverage = sample_coverage,
    estimator = estimator,
    as_numeric = FALSE,
    check_arguments = FALSE
  )

  return(
    # Make a tibble with site, estimator and sample-coverage
    tibble::tibble(
      # Restore non-species columns
      x[colnames(x) %in% non_species_columns],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, size_list)
    )
  )
}


#' Departure of actual sample coverage from target coverage
#'
#' Helper for [coverage_2_size].
#'
#' @param size The size of the sample. Adjusted to minimize `delta()`.
#' @param target_coverage The sample coverage to reach by adjusting size.
#'
#' @returns The departure of actual sample coverage from target coverage.
#' @noRd
#'
chao_delta <- function(
    abd,
    size,
    target_coverage) {
  abs(
    coverage.numeric(
      x = abd,
      estimator = "Chao",
      level = size,
      as_numeric = TRUE,
      check_arguments = FALSE
    ) - target_coverage
  )
}
