#' Abundance Frequency Count of a Community
#'
#' Count the number of species observed the same number of times.
#'
#' The Abundance Frequency Count \insertCite{Chao2015}{divent} is the number
#' of species observed each number of times.
#' It is a way to summarize the species distribution.
#'
#' It can be estimated at a specified level of interpolation or extrapolation.
#' Extrapolation relies on the estimation of the estimation of the asymptotic
#' distribution of the community by [probabilities] and eq. (5)
#' of \insertCite{Chao2014}{divent}.
#'
#' @inheritParams check_divent_args
#' @param abd A numeric vector containing species abundances.
#' @param richness_estimator A string containing an estimator recognized by
#' [div_richness] to evaluate the total number of species in [probabilities].
#' Used only for extrapolation.
#'
#' @returns A two-column tibble. The first column contains the number of observations,
#' the second one the number of species observed this number of times.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' abd_freq_count(as.numeric(paracou_6_abd[1, ]))
abd_freq_count <- function(
    abd,
    level = NULL,
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
    jack_alpha = 0.05,
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    check_arguments = TRUE) {

  # Check arguments
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  richness_estimator <- match.arg(richness_estimator)
  if (any(check_arguments)) check_divent_args()

  # Convert to integer values
  if (length(abd) == 0) {
    # abd may contain no species. Return NA.
    return(NA)
  } else {
    abd_int <- round(abd)
  }
  if (any(abs(abd_int - abd) > 2 * .Machine$double.eps)) {
    cli::cli_alert_warning(
      "The abundance frequency count requires integer abundances."
    )
    cli::cli_alert(
      "Abundances have been rounded."
    )
  }

  # Eliminate 0
  abd_int <- abd_int[abd_int > 0]
  # Calculate abundance distribution
  abd_distribution <- tapply(abd_int, INDEX = abd_int, FUN = length)

  if (is.null(level)) {
    # No extrapolation. Prepare a two-column matrix ----
    return(
      tibble::tibble(
        abundance = as.numeric(names(abd_distribution)),
        number_of_species = abd_distribution
      )
    )
  } else {
    # If level is coverage, get size
    if (level < 1) {
      level <- coverage_to_size.numeric(
        abd,
        sample_coverage = level,
        estimator = coverage_estimator,
        as_numeric  = TRUE,
        check_arguments = FALSE
      )
    }
    sample_size <- sum(abd_int)
    if (level <= sample_size) {
      # Interpolation ----
      s_nu <- vapply(
        seq_len(level),
        function(nu) {
          sum(
            exp(
              lchoose(abd_int, nu) +
              lchoose(sample_size - abd_int, level - nu) -
              lchoose(sample_size, level)
            )
          )
        },
        FUN.VALUE = 0
      )
    } else {
      # Extrapolation ----
      if (length(abd_int) == 1) {
        # Single species: general formula won't work: log(1-prob_s_0)
        s_nu <- c(rep(0, level - 1), 1)
      } else {
        # Unveil the full distribution
        prob_s_0 <- probabilities.numeric(
          abd_int,
          estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator,
          jack_alpha = jack_alpha,
          jack_max = jack_max,
          check_arguments = TRUE
        )
        # Extrapolate
        s_nu <- vapply(
          seq_len(level),
          function(nu) {
            sum(
              exp(
                lchoose(level, nu) +
                nu * log(prob_s_0) +
                (level - nu) * log(1 - prob_s_0)
              )
            )
          },
          FUN.VALUE = 0
        )
      }
    }
    # Return a tibble with all possible abundances
    return(
      tibble::tibble(
        abundance = seq_len(level),
        number_of_species = s_nu)
    )
  }
}
