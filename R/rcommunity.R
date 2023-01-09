#' Random communities
#' 
#' Draw random communities according to a probability distribution.
#' 
#' Communities of fixed `size` are drawn in a multinomial distribution according 
#' to the distribution of probabilities provided by `prob`.
#' An abundance vector `abd` may be used instead of probabilities, 
#' then `size` is by default the total number of individuals in the vector. 
#' Random communities can be built by drawing in a multinomial law following 
#' \insertCite{Marcon2012a;textual}{divent}, or trying to estimate the 
#' distribution of the actual community with [as_probabilities]. 
#' If `bootstrap` is "Chao2013", the distribution is estimated by a single 
#' parameter model and unobserved species are given equal probabilities. 
#' If `bootstrap` is "Chao2015", a two-parameter model is used and unobserved 
#' species follow a geometric distribution.
#' 
#' Alternatively, the probabilities may be drawn following a classical 
#' distribution: either lognormal ("lnorm")  \insertCite{Preston1948}{divent} 
#' with given standard deviation (`sd_lnorm`; note that the mean is actually 
#' a normalizing constant. Its values is set equal to 0 for the simulation of 
#' the normal distribution of unnormalized log-abundances), log-series ("lseries")
#' \insertCite{Fisher1943}{divent} with parameter `alpha_lseries`, geometric 
#' ("geom") one \insertCite{Motomura1932}{divent} with parameter `prob_geom`, 
#' or broken stick ("bstick") \insertCite{MacArthur1957}{divent}. 
#' The number of simulated species is fixed by `species_number`, except for 
#' "lseries" where it is obtained from `alpha_lseries` and `size`: 
#' \eqn{S=\alpha \ln(1 + size / \alpha)}.
#' 
#' Log-normal, log-series and broken-stick distributions are stochastic. 
#' The geometric distribution is completely determined by its parameters.
#'
#' @inheritParams check_divent_args
#' @param n The number of communities to draw.
#' @param size The number of individuals to draw in each community.
#' @param prob A numeric vector containing probabilities.
#' @param abd A numeric vector containing abundances.
#' @param bootstrap The method used to obtain the probabilities to generate 
#' bootstrapped communities from observed abundances. 
#' If "Marcon2012", the probabilities are simply the abundances divided by the total
#' number of individuals \insertCite{Marcon2012a}{divent}. 
#' If "Chao2013" or "Chao2015" (by default), a more sophisticated approach is used 
#' (see [as_probabilities]) following\insertCite{Chao2013;textual}{divent} or 
#' \insertCite{Chao2015;textual}{divent}.
#' @param species_number The number of species.
#' @param distribution The distribution of species abundances.
#' May be "lnorm" (log-normal), "lseries" (log-series), "geom" (geometric) or 
#' "bstick" (broken stick).
#' @param sd_lnorm The simulated log-normal distribution standard deviation. 
#' This is the standard deviation on the log scale.
#' @param prob_geom The proportion of ressources taken by successive species 
#' of the geometric distribution.
#' @param alpha_lseries Fisher's alpha in the log-series distribution.
#'
#' @return An object of class [abundances].
#' @export
#'
#' @examples
#' # Generate a community made of 100000 individuals among 300 species and fit it
#' abundances <- rcommunity(n = 1, size = 1E5, 
#'   species_number = 300, distribution = "lnorm")
#' autoplot(abundances)
rcommunity <- function(
    n,
    size = sum(abd),
    prob = NULL,
    abd = NULL,
    bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
    species_number = 300,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    sd_lnorm = 1,
    prob_geom = 0.1,
    alpha_lseries = 40,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (!is.null(prob) & !is.null(abd)) {
      stop("'prob' and 'abd' can't be both given.")
    }
  }
  bootstrap <- match.arg(bootstrap) 
  distribution <- match.arg(distribution) 
  coverage_estimator <- match.arg(coverage_estimator)
  
  the_prob <- the_abd <- NULL
  
  # Generate a distribution (prob and abd are null) or use prob or build probabilities from abd
  if (is.null(prob) & is.null(abd)) {
    # Draw in a distribution
    name <- distribution
    if (distribution == "lseries") {
      # Draw probabilities, except for logseries: draw abundances.
      the_abd <- t(
        replicate(
          # Number of communities
          n,
          # Draw each species separately
          replicate(
            # Number of species (Fisher's formula)
            round(-alpha_lseries * log(alpha_lseries / (size + alpha_lseries))),
            # Internal function to draw the number of individuals of a species
            abd_lseries(size, alpha_lseries)
          )
        )
      )
    } else {
      # Other distributions: draw probabilities
      the_prob <- switch(distribution,
        geom = prob_geom / (1 - (1 - prob_geom) ^ species_number) * (1 - prob_geom) ^ (0:(species_number - 1)),
        lnorm = (n_lnorm <- stats::rlnorm(species_number, 0, sd_lnorm)) / sum(n_lnorm),
        bstick = c(cuts <- sort(stats::runif(species_number - 1)), 1) - c(0, cuts)
      )
    }
  } else {
    # Simulation from a distribution
    name <- ifelse(is.null(prob), "prob", "abd")
  }
  if (!is.null(prob)) {
    # Probabilities are given
    the_prob <- prob[prob != 0]
  }
  if (!is.null(abd)) {
    # Subsample in given abundances. Generate probabilities according to the chosen method.
    if (bootstrap == "Chao2015") {
      the_prob <- probabilities.numeric(
        abd,
        estimator = "Chao2015",
        unveiling = "geometric",
        coverage_estimator = coverage_estimator,
        as_numeric = TRUE,
        check_arguments = FALSE)
    }
    if (bootstrap == "Chao2013") {
      the_prob <- probabilities.numeric(
        abd, 
        estimator = "Chao2013", 
        unveiling = "uniform", 
        coverage_estimator = coverage_estimator,
        as_numeric = TRUE,
        check_arguments = FALSE)
    }
    if (bootstrap == "Marcon2012") {
      the_prob <- abd / sum(abd)
    }
  }
  
  # Generate communities according to probabilities
  if (is.null(the_abd)) {
    # Draw multinomial samples from probabilities except if abundances have already been obtained (e.g.: lseries)
    the_abd <- t(stats::rmultinom(n, size, the_prob))
  }

  return(
    as_abundances.matrix(
      the_abd,
      names = paste(
        name, 
        1:n,
        sep = "_"
      ),
      round = FALSE,
      check_arguments = FALSE
    )
  )
}


#' Abundance in a log-series
#' 
#' Abundance of a species in a logseries distribution of given 
#' size and Fisher's alpha.
#' 
#' Adapted from Dan Lunn, http://www.stats.ox.ac.uk/~dlunn/BS1_05/BS1_Rcode.pdf
#' 
#' @param size The number of individuals in the community.
#' @param alpha_lseries The value of Fisher's alpha.
#'
#' @return The number of individuals of the species.
#' @noRd
#' 
abd_lseries <- function(size, alpha_lseries) {
  # Fisher's x is log-series 1-theta
  x <- size / (size + alpha_lseries)
  # Draw a random number between 0 and 1
  u <- stats::runif(1)
  # k is the number of individuals to draw
  k <- 1
  # Calculate the probability at k=1
  p <- -x / log(1 - x)
  # Store it in the distribution function
  p_cumulated <- p
  # Repeat while the cumulated probability is below u
  while (p_cumulated <= u) {
    # Probability at k+1 obtained from that at k
    p <- p * k * x / (k + 1)
    # Increment k
    k <- k + 1
    # Increment the cumulated probability
    p_cumulated <- p_cumulated + p
  }
  return(k)
}
