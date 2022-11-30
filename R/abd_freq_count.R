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
#' @param abd A numeric vector containing species abundances.
#' @param level The level of interpolation or extrapolation. 
#' It may be an an arbitrary sample size (an integer) or a sample coverage
#' (a number between 0 and 1).
#' @param probability_estimator A string containing one of the possible estimators
#' of the probability distribution (see [probabilities]). 
#' Used only for extrapolation.
#' @param unveiling A string containing one of the possible unveiling methods 
#' to estimate the probabilities of the unobserved species (see [probabilities]).
#' Used only for extrapolation.
#' @param richness_estimator A string containing an estimator recognized by 
#' [div_richness] to evaluate the total number of species in [probabilities]. 
#' Used only for extrapolation.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#'
#' @return A two-column tibble. The first column contains the number of observations, 
#' the second one the number of species observed this number of times.
#' @export
#'
#' @examples
#' abd_freq_count(paracou_6_abd[1, -(1:2)])
abd_freq_count <- function (
    abd,
    level = NULL, 
    probability_estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy"),
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator)

  # Convert to integers
  if (length(abd) == 0) {
    # abd may contain no species. Return NA.
    return(NA)
  } else {
    abd_int <- as.integer(round(abd))    
  }
  if (any(abs(abd_int - abd)> 2 * .Machine$double.eps))
    warning ("The abundance frequency count requires integer abundances. Abundances have been rounded.")
  
  # Eliminate 0
  abd_int <- abd_int[abd_int > 0]
  # Calculate abundance distribution
  abd_distribution <- tapply(abd_int, abd_int, length)
  
  if (is.null(level)) {
    # No extrapolation. Prepare a two-column matrix ----
    return(
      tibble::tibble(
        abundance = as.integer(names(abd_distribution)), 
        number_of_species = abd_distribution)
    )
  } else {
    # If level is coverage, get size
    if (level < 1) {
      level <- coverage_to_size.numeric(
        abd, 
        sample_coverage = level,
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
        FUN.VALUE=0.0
      )
    } else {
      # Extrapolation ----
      if (length(abd_int) == 1) {
        # Single species: general formula won't work: log(1-prob_s_0)
        s_nu <- c(rep(0, level - 1), 1)
      } else {
        # Unveil the full distribution
        prob_s_0 <- probabilities(
          abd_int,
          estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator, 
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
          } , 
          FUN.VALUE=0.0
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
