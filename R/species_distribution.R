#' Species Distributions
#' 
#' A Species Distribution is a [tibble::tibble] containing species abundances or probabilities.
#' 
#' `species_distribution` objects include `abundances` and `probabilities` objects.
#' 
#' `as_species_distribution()`, `as_abundances()`  and `as_probabilities` format 
#' the numeric, matrix or dataframe `x` so that appropriate 
#' versions of community functions (generic methods such as [plot] or 
#' [div_richness]) are applied. 
#' Abundance values are rounded (by default) to the nearest integer.
#' 
#' `as_probabilities()` normalizes the vector `x` so that it sums to 1. It gives
#' the same output as `probabilities()` with `estimator = "naive"`.
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
#' TODO: These functions can be applied to data frames to calculate the joint diversity \insertCite{Gregorius2010}{divent}.
#'  
#' TODO: `species_distribution` objects objects can be plotted by `plot` and `autoplot()`.
#'
#' @param x An object.
#' @param ... Additional arguments to be passed to [plot]. Unused elsewhere.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#' 
#' @examples
#' # Paracou data
#' is_species_distribution(paracou_6_abd)
#' # Whittaker plot fitted by a log-normal distribution TODO
#' # autoplot(asm_paracou_6, distribution = "lognormal")
#' @references
#' \insertAllCited{}
#' 
#' @import tibble
#' 
#' @name species_distribution
NULL


#  Species Distribution ----

#' @rdname species_distribution
#'
#' @param names The names of the species distributions.
#' @param weights The weights of the sites of the species distributions.
#' 
#' @export
species_distribution <- function(
    x, 
    names = NULL, 
    weights = NULL, 
    check_arguments = TRUE) {

  # Check the data ----
  if (check_arguments) check_divent_args()
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (length(dim(x)) > 2) stop("'x' may be a vector or a matrix")
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  
  # Build a tibble from the data ----
  if (is.vector(x)) {
    ## Single distribution ----
    if (is.null(names(x))) {
      ### Columns: add default species names such as sp_1 ----
      names(x) <- paste(
        "sp", 
        formatC(seq_along(x), width = ceiling(log10(length(x))), flag = "0"),
        sep = "_"
      )
    }
    if (length(names) != 1) {
      ### Rows: Add a site name ----
      names <- paste("site", round(stats::runif(1)*.Machine$integer.max), sep="_")
    }
    # Build a tibble
    distribution <- tibble::as_tibble_row(c(site = names, x))
    
  } else {
    ## Several distributions ----
    if (is.null(colnames(x))) {
      ### Columns: add default species names such as sp_1 ----
      colnames(x) <- paste(
        "sp", 
        formatC(seq_len(ncol(x)), width = ceiling(log10(ncol(x))), flag = "0"),
        sep = "_"
      )
    }
    # Build a tibble
    distribution <- tibble::as_tibble(x, rownames = "site")
    ### Rows: site names = names or matrix row names or default ----
    if (!is.null(names)) {
      # site = names if the size matches
      if (length(names) == nrow(x)) {
        distribution$site <- names
      } else {
        stop("The length of 'names' must match the number of lines of the data matrix.")
      }
    } else {
      # names is null...
      if (is.null(row.names(x))) {
        # ...and no row names: set default names such as site_1
        distribution$site <- paste(
          "site", 
          formatC(seq_len(nrow(x)), width = ceiling(log10(nrow(x))), flag = "0"),
          sep = "_"
        )
      }
    }
    ### Rows: site weights ----
    if (!is.null(weights)) {
      # site = weights if the size matches
      if (length(weights) == nrow(x)) {
        distribution <- tibble::add_column(
          distribution, 
          weight = rowSums(x),
          .after = "site"
        )
      } else {
        stop("The length of 'weights' must match the number of lines of the data matrix.")
      }
    } else {
      # Weights are the number of individuals
      distribution <- tibble::add_column(
        distribution, 
        weight = rowSums(x),
        .after = "site"
      )
    }
  }
  
  # Set the class and return ----
  class(distribution) <- c("species_distribution", class(distribution))
  return(distribution)
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution <- function(x, ...) {
  UseMethod("as_species_distribution")
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.numeric <- function(
    x,
    ...,
    check_arguments = TRUE) {
  
  return(
    species_distribution(
      x, 
      ..., 
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.matrix <- function(
    x, 
    names = NULL, 
    weights = NULL, 
    ...,
    check_arguments = TRUE) {
  
  return(
    species_distribution(
      x, 
      names = names, 
      weights = weights, 
      ..., 
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  # Check the data
  if (any(x < 0)) stop("All numeric values of the dataframe must be positive.")
  
  # Build a tibble
  distribution <- tibble::as_tibble(x)
  
  # The first column should be "site"
  if (!"site" %in% colnames(distribution)) {
    distribution <- tibble::add_column(
      distribution, 
      site = paste(
        "site", 
        formatC(
          seq_len(nrow(distribution)), 
          width = ceiling(log10(nrow(distribution))), 
          flag = "0"
        ),
        sep = "_"
      ),
      .before = 1
    )
  }
  
  # The second column should be "weight"
  if (!"weight" %in% colnames(distribution)) {
    distribution <- tibble::add_column(
      distribution, 
      weight = rowSums(distribution[colnames(distribution) != "site"]),
      .after = "site"
    )
  }

  # Set the class and return
  class(distribution) <- c("species_distribution", class(distribution))
  return(distribution)
}


#' @rdname species_distribution
#'
#' @export
is_species_distribution <- function(x) {
  inherits(x, "species_distribution")
}

#  Probabilities ----

#' @rdname species_distribution
#'
#' @export
probabilities <- function(x, ...) {
  UseMethod("probabilities")
}


#' @rdname species_distribution
#'
#' @param abd A numeric vector containing abundances. It may be named to track
#' species names.
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
#' @param jack_max The highest jackknife order allowed. Default is 10. 
#' Allowed values are between 1 and 10.
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
    jack_max = 10, 
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q = 0,
    ...,
    check_arguments = TRUE) {
  
  # Check the data ----
  if (check_arguments) check_divent_args()
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
    prob <- x / sum(x)
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
    sample_coverage <- coverage(abd, estimator = coverage_estimator)$coverage
    if (
      estimator == "Chao2015" | 
      unveiling == "Chao2015" | 
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
      s_0 <- div_richness.numeric(abd, estimator = "jackknife", check_arguments = FALSE)$richness - s_obs
      # Estimate the number of unobserved species by iterations
# TODO : activate ent_tsallis
#      ent_target <- ent_tsallis(abd, q = q, estimator = "naive", check_arguments = FALSE)$entropy
      ent_target <- div_richness(abd, estimator = "naive", check_arguments = FALSE)$richness
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
        div_richness(
          abd, 
          estimator = richness_estimator, 
          jack_max = jack_max,
          check_arguments = FALSE
        )$richness
      )
      s_0 <- s_est - s_obs
    }
    
    ## Distribution of unobserved species ----
    if (s_0) {
      if (unveiling == "none") {
        prob <- prob_tuned
      } else {
        prob <- c(
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
      prob <- prob_tuned
    }
  }
  
  # Set the class and return ----
  probabilities <- as_species_distribution(prob)
  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#'
#' @export
probabilities.abundances <- function(
    x, 
    estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"), 
    jack_max = 10, 
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q = 0,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
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
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    q = q,
    check_arguments = FALSE
  )
  
  # Bind the rows
  probabilities <- dplyr::bind_rows(probabilities_list)
  # Restore the site names
  if (!is.null(x$site)) probabilities$site <- x$site
  # Replace NA's due to binding by zeros
  probabilities <- dplyr::mutate_all(probabilities,  ~replace(., is.na(.), 0))
  return(probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities <- function(x, ...) {
  UseMethod("as_probabilities")
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.numeric <- function(
    x,
    ..., 
    check_arguments = TRUE) {

  if (any(x < 0)) stop("Species probabilities must be positive.")

  # Normalize to 1
  prob <- x / sum(x)
  probabilities <- as_species_distribution(
    prob, 
    ...,
    check_arguments = check_arguments
  )

  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#' 
#' @export
as_probabilities.matrix <- function(
    x,
    names = NULL, 
    weights = NULL, 
    ...,
    check_arguments = TRUE) {
  
  probabilities <- as_species_distribution.matrix(
    x, 
    names = names, 
    weights = weights,
    ..., 
    check_arguments = check_arguments
  )
  
  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  probabilities <- as_species_distribution.data.frame(
    x, 
    ..., 
    check_arguments = check_arguments
  )
  
  class(probabilities) <- c("probabilities", class(probabilities))
  return(probabilities)
}


#' @rdname species_distribution
#' 
#' @export
is_probabilities <- function(x) {
  inherits(x, "probabilities")
}


#  Abundances ----

#' @rdname species_distribution
#'
#' @param round If `TRUE`, the values of `x` are converted to integers.
#' 
#' @export
abundances <- function(
    x,
    round = TRUE,
    names = NULL, 
    weights = NULL, 
    check_arguments = TRUE) {
  
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (any(x < 0)) stop("Species abundances must be positive.")
  if (round) {
    # Add 0.5 before changing mode to round rather than taking the floor
    x <- x + 0.5
    mode(x) <- "integer"
  }
  
  abundances <- species_distribution(
    x,     
    names = names, 
    weights = weights, 
    check_arguments = check_arguments
  )

  class(abundances) <- c("abundances", class(abundances))
  return(abundances)    

}

#' @rdname species_distribution
#' 
#' @export
as_abundances <- function(x, ...) {
  UseMethod("as_abundances")
}


#' @rdname species_distribution
#' 
#' @export
as_abundances.numeric <- function(
    x,
    round = TRUE, 
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  if (any(x < 0)) stop("Species abundances must be positive.")
  
  if (round) {
    x <- as.integer(round(x))
  }

  abundances <- as_species_distribution(x, ..., check_arguments = FALSE)

  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#' 
#' @export
as_abundances.matrix <- function(
    x,
    round = TRUE,
    names = NULL, 
    weights = NULL, 
    ...,
    check_arguments = TRUE) {
  
  if (round) {
    # Add 0.5 before changing mode to round rather than taking the floor
    x <- x + 0.5
    mode(x) <- "integer"
  }
  
  abundances <- as_species_distribution.matrix(
    x, 
    names = names, 
    weights = weights,
    ..., 
    check_arguments = check_arguments
  )
  
  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.data.frame <- function(
    x, 
    ...,
    check_arguments = TRUE) {
  
  abundances <- as_species_distribution.data.frame(
    x, 
    ..., 
    check_arguments = check_arguments
  )
  
  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#' 
#' @export
is_abundances <- function(x) {
  inherits(x, "abundances")
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

  prob_s_0 <- NA
  if (unveiling == "geometric") {
    if (s_0 == 1) {
      # A single unobserved species
      prob_s_0 <- 1 - sample_coverage
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
      prob_s_0 <- alpha * beta^i
      # Sometimes fails when the distribution is very uneven (sometimes r < 1) 
      # Then, fall back to the uniform distribution
      if (any(is.na(prob_s_0)) | any(prob_s_0 <= 0)) {
        unveiling <- "uniform"
      }
    }
  }      
  if (unveiling == "uniform") {
    # Add s_0 unobserved species with equal probabilities
    prob_s_0 <- rep((1 - sum(prob_tuned)) / s_0, s_0)
  }
  if (any(is.na(prob_s_0))) {
    warning("Unveiling method was not recognized")
    return(NA)
  } else {
    names(prob_s_0) <- paste("UnobsSp", seq_along(prob_s_0), sep = "")
    return(prob_s_0)
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
    ent_bias <- abs(
      sum(
        -seq_len(sample_size) / sample_size * 
          log(seq_len(sample_size) / sample_size) * s_nu
      ) 
      - ent_target
    )
  } else {
    ent_bias <- abs(
      (sum((seq_len(sample_size)/sample_size)^q * s_nu) - 1) / (1 - q) 
      - ent_target
    )
  }
  return(ent_bias)
}
