#' Species Distributions
#' 
#' A Species Distribution is a [tibble::tibble] containing species abundances or probabilities.
#' 
#' `species_distribution` objects include `abundances` and `probabilities` objects.
#' 
#' `as_abundances()` formats the numeric or integer `x` so that appropriate versions of community functions (generic methods such as Diversity) are applied. Abundance values are rounded (by default) to the nearest integer.
#' 
#' `as_probabilities()` normalizes the vector `x` so that it sums to 1. 
#' If the `estimator` is not "naive", the observed abundance distribution is used 
#' to estimate the actual species distribution. The list of species will be changed:
#' zero-abundance species will be cleared, and some unobserved species will be added. 
#' First, observed species probabilities are estimated folllowing 
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
#' `species_distribution` objects objects can be plotted by `plot` and `autoplot()`.
#'
#' @param x An object.
#' @param name The name of the species distribution.
#' @param ... Additional arguments to be passed to [plot]. Unused elsewhere.
#' @param check_arguments If `TRUE`, check the arguments of the function. 
#' 
#' @examples
#' # Paracou data TODO
#' is_species_distribution(asm_paracou_6)
#' # Whittaker plot fitted by a log-normal distribution
#' autoplot(asm_paracou_6, distribution = "lognormal")
#' @references
#' \insertAllCited{}
#' 
#' @import tibble
#' 
#' @name species_distribution
NULL


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
    name = "", 
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  # Check the data
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  
  # Add default names such as sp_1
  if (is.null(names(x))) {
    names(x) <- paste(
      "sp", 
      formatC(1:length(x), width = ceiling(log10(length(x))), flag = "0"),
      sep = "_"
    )
  }
  
  # Build a tibble
  distribution <- tibble::as_tibble_row(c(site = name, x))
  
  # Set the class and return
  class(distribution) <- c("species_distribution", class(distribution))
  return(distribution)
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.integer <- function(
    x, 
    name = "", 
    ...,
    check_arguments = TRUE) {
  return(
    as_species_distribution.numeric(
      x, 
      name = name, 
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
    ...,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  # Check the data
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  
  # Columns: add default names such as sp_1
  if (is.null(colnames(x))) {
    colnames(x) <- paste(
      "sp", 
      formatC(1:ncol(x), width = ceiling(log10(ncol(x))), flag = "0"),
      sep = "_"
    )
  }
  
  # Build a tibble
  distribution <- tibble::as_tibble(x, rownames = "site")
  # Weights are the number of individuals
  distribution <- tibble::add_column(
    distribution, 
    weight = rowSums(x),
    .after = "site"
  )
  
  # Lines: add default names such as site_1
  if (!is.null(names)) {
    # site = names if the size matches
    if (length(names) == nrow(x)) {
      distribution$site <- names
    } else {
      stop("The length of names must match the number of lines of the data matrix.")
    }
  } else {
    # names is null...
    if (is.null(row.names(x))) {
      # ...and no row names: set default names
      distribution$site <- paste(
        "site", 
        formatC(1:nrow(x), width = ceiling(log10(nrow(x))), flag = "0"),
        sep = "_"
      )
    }
  }
  
  # Set the class and return
  class(distribution) <- c("species_distribution", class(distribution))
  return(distribution)
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
          1:nrow(distribution), 
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


#' @rdname species_distribution
#'
#' @export
as_probabilities <- function(x, ...) {
  UseMethod("as_probabilities")
}


#' @rdname species_distribution
#'
#' @param x 
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
as_probabilities.numeric <- function(
    x, 
    estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "rarefy"), 
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
  
  
  if (estimator == "naive") {
    # Just normalize so that x sums to 1
    x <- x / sum(x)
  } else {
    # Integer abundances are required by all non-naive estimators
    if (!is_integer_values(x)) {
      warning(
        "Integer abundance values are required to estimate community probabilities. Abundances have been rounded."
      )
    }
    x_int <- round(x)
    
    # Eliminate 0 and calculate elementary statistics
    abundances <- x_int[x_int > 0]
    species_number <- length(abundances)
    sample_size <- sum(abundances)
    probabilities <- abundances / sample_size
    # Sample coverage
    sample_coverage <- coverage(abundances, estimator = coverage_estimator)
    if (
        estimator == "Chao2015" | 
        unveiling == "Chao2015" | 
        richness_estimator == "Rarefy") {
      # Sample coverage of order 2 required
      singletons <- sum(abundances == 1)
      doubletons <- sum(abundances == 2)
      if (doubletons == 0) {
        singletons <- max(singletons - 1, 0)
        doubletons <- 1
      }
      tripletons <- max(sum(abundances == 3), 1)
      # 1 minus sample coverage (i.e. Coverage Deficit) of order 2
      coverage_deficit_2 <-
        doubletons / choose(sample_size, 2) * 
        (
          (sample_size - 2) * doubletons / 
            ((sample_size - 2) * doubletons + 3 * tripletons)
        )^2
    }
    
    # Tune the probabilities of observed species
    if (sample_coverage == 0 | sample_coverage == 1) {
      # Sample coverage equal to 1, do not tune. If 0, unable to tune.
      probabilities_tuned <- probabilities
    } else {
      probabilities_tuned <- NA
      if (estimator == "ChaoShen") {
        probabilities_tuned <- sample_coverage * probabilities
      }
      if (estimator == "Chao2013") {
        # Single parameter estimation, Chao et al. (2013)
        denominator <- sum(probabilities * (1 - probabilities)^sample_size)
        if (denominator == 0) {
          # N is too big so denominator equals 0. Just multiply by C.
          probabilities_tuned <- sample_coverage * probabilities
        } else {
          # General case
          lambda <- (1 - sample_coverage) / denominator
          probabilities_tuned <- probabilities * (1 - lambda * (1 -probabilities)^sample_size)
        }      
      } 
      if (estimator == "Chao2015")  {
        # Two parameters, Chao et al. (2015). 
        # Estimate theta. Set it to 1 if impossible
        theta <- tryCatch(
          stats::optimize(
            theta_solve, 
            interval=c(0,1), 
            probabilities, 
            abundances, 
            sample_size, 
            sample_coverage, 
            coverage_deficit_2
          )$min, 
          error = function(e) {1}
        )
        lambda <- (1 - sample_coverage) / sum(probabilities * exp(-theta * abundances))
        probabilities_tuned <- probabilities * (1 - lambda * exp(-theta * abundances))
      }
    }
    names(PsTuned) <- names(spD[spD > 0])
    
    # Estimate the number of unobserved species
    if (RCorrection == "Rarefy") {
      if (Unveiling == "none")
        stop("Arguments RCorrection='Rarefy' and Unveiling='None' are not compatible")
      # Estimation of the number of unobserved species to initialize optimization
      S0 <- bcRichness(Ns, Correction="Jackknife") - S
      # Estimate the number of unobserved species by iterations
      Target <- Tsallis(Ns, q=q, Correction="None", check_arguments = FALSE)
      S0 <- round(tryCatch(stats::optimize(rarefaction_bias, interval=c(0, 2*S0), Ns, PsTuned, C, CD2, q, Unveiling, Target)$minimum,
                           error = function(e) {S0}))
    } else {
      Sestimate <- ceiling(bcRichness(Ns, Correction=RCorrection, JackOver=JackOver, JackMax=JackMax))
      S0 <- Sestimate - S
    }
    
    # Distribution of unobserved species
    if (S0) {
      if (Unveiling == "None") {
        spD <- PsTuned
      } else {
        spD <- c(PsTuned, estimate_Ps0(Unveiling, PsTuned, S0, C, CD2))
      }
    } else {
      spD <- PsTuned
    }
    spD <- as_species_distribution(spD, ...)
  }
  class(spD) <- c("probabilities", class(spD))
  return(spD)
}


#' @rdname species_distribution
#' 
#' @export
as_probabilities.integer <- function(
    x, 
    estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "rarefy"), 
    jack_max = 10, 
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    q = 0, 
    ..., 
    check_arguments = TRUE) {
  
  return(
    as_probabilities.numeric(
      x, 
      estimator = estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator, 
      jack_max = jack_max, 
      coverage_estimator = coverage_estimator,
      q = q, 
      ..., 
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#' 
#' @export
is_probabilities <- function(x) {
  inherits(x, "probabilities")
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
    name = "",
    round = TRUE, 
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()

  if (round) {
    x <- as.integer(round(x))
  }

  abundances <- as_species_distribution(x, name = name, ..., check_arguments = FALSE)

  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}


#' @rdname species_distribution
#' 
#' @export
as_abundances.integer <- function(
    x, 
    name = "", 
    ..., 
    check_arguments = TRUE) {

  return(
    as_abundances.numeric(
      x, 
      name = name, 
      round = FALSE, 
      ...,
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#' 
#' @export
as_abundances.matrix <- function(
    x,
    names = "",
    round = TRUE, 
    ...,
    check_arguments = TRUE) {
  
  if (round) {
    # Add 0.5 before changing mode to round rather than taking the floor
    x <- x + 0.5
    mode(x) <- "integer"
  x}
  
  abundances <- as_species_distribution.matrix(
    x, 
    names = names, 
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


# Utilities for as_probabilities.numeric. Not exported. ####
# Solve the theta parameter of Chao et al. (2015)
theta_solve <- function(theta, Ps, Ns, N, C, CD2){
  # Code inspired from JADE function DetAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
  lambda <- (1-C) / sum(Ps * exp(-theta*Ns))
  return(abs(sum((Ps * (1 - lambda * exp(-theta*Ns)))^2) - sum(choose(Ns,2)/choose(N,2)) + CD2))
}
# Solve the beta parameter of Chao et al. (2015)
beta_solve <- function(beta, r, i){
  # Code inspired from JADE function UndAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
  return(abs(sum(beta^i)^2 / sum((beta^i)^2) - r))
}
# Unobserved species distribution
estimate_Ps0 <- function(Unveiling, PsTuned, S0, C, CD2){
  Ps0 <- NA
  if (Unveiling == "geom") {
    if (S0 == 1) {
      # A single unobserved species
      Ps0 <- 1-C
    } else {
      r <- (1-C)^2/CD2
      i <- seq_len(S0)
      beta <-  tryCatch(stats::optimize(beta_solve, lower=(r-1)/(r+1), upper=1, tol=.Machine$double.eps, r, i)$min, 
                        error = function(e) {(r-1)/(r+1)})
      alpha <- (1-C) / sum(beta^i)
      Ps0 <- alpha * beta^i
      # Sometimes fails when the distribution is very uneven (sometimes r < 1) 
      # Then, go back to the uniform distribution
      if (any(is_na(Ps0)) | any(Ps0 <= 0)) Unveiling <- "unif"
    }
  }      
  if (Unveiling == "unif") {
    # Add S0 unobserved species with equal probabilities
    Ps0 <- rep((1-sum(PsTuned))/S0, S0)
  }
  if (any(is_na(Ps0))) {
    warning("Unveiling method was not recognized")
    return(NA)
  } else {
    names(Ps0) <- paste("UnobsSp", seq_along(length(Ps0)), sep="")
    return(Ps0)
  }         
}
# Rarefaction bias
rarefaction_bias <- function(S0, Ns, PsTuned, C, CD2, q, Unveiling, Target) {
  Ns <- Ns[Ns>0]
  N <- sum(Ns)
  # Unobserved species
  Ps0 <- estimate_Ps0(Unveiling, PsTuned, S0, C, CD2)
  # Full distribution of probabilities
  Ps <- c(PsTuned, Ps0)
  # AbdFreqCount at Level = N
  Sn <- vapply(seq_len(N), function(nu) sum(exp(lchoose(N, nu) + nu*log(Ps) + (N-nu)*log(1-Ps))), FUN.VALUE=0.0)
  # Get Entropy at Level=N and calculate the bias
  if (q == 1) {
    Bias <- abs(sum(-seq_len(N)/N * log(seq_len(N)/N) * Sn) - Target)
  } else {
    Bias <- abs((sum((seq_len(N)/N)^q * Sn) - 1) / (1-q) - Target)
  }
  return(Bias)
}
# end of utilities ####