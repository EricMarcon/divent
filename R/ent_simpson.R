#' Simpson's Entropy of a Community
#' 
#' Estimate the entropy \insertCite{Simpson1949}{divent} of species from abundance
#' or probability data.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#' See [div_hill] for non-specific estimators.
#' 
#' Simpson-specific estimator is from \insertCite{Lande1996;textual}{divent}.
#' 
#' Entropy can be estimated at a specified level of interpolation or 
#' extrapolation, either a chosen sample size or sample coverage 
#' \insertCite{Chao2014}{divent}, rather than its asymptotic value.
#' See [ent_accum] for details.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @return A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # Entropy of each community
#' ent_simpson(paracou_6_abd)
#' # gamma entropy
#' ent_simpson(paracou_6_abd, gamma = TRUE)
#' 
#' # At 80% coverage
#' ent_simpson(paracou_6_abd, level = 0.8)
#' 
#' @name ent_simpson
NULL


#' @rdname ent_simpson
#'
#' @export
ent_simpson <- function(x, ...) {
  UseMethod("ent_simpson")
}


#' @rdname ent_simpson
#'
#' @param estimator An estimator of entropy.
#' 
#' @export
ent_simpson.numeric <- function(
    x, 
    estimator = c("Lande",
                  "UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    if (!is.null(level)) stop("Entropy can't be estimated at a level without the abundance of species.")
    # Probabilities sum to 1, allowing rounding error
    the_entropy <- 1 - sum(x^2)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive", 
          order = 2,
          entropy = the_entropy
        )
      )  
    }
  }
  
  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)

  # Asymptotic estimator ----
  if (is.null(level)) {
    ## Exit if x contains no or a single species ----
    if (length(abd) < 2) {
      if (length(abd) == 0) {
        if (as_numeric) {
          return(NA)
        } else {
          return(
            tibble::tibble_row(
              estimator = "No Species", 
              order = 2,
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
              order = 2,
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
    
    if (estimator == "Lande") {
      # Lande's estimator ----
      the_entropy <- 1 - sum(abd * (abd - 1) / sample_size / (sample_size - 1))
      if (as_numeric) {
        return(the_entropy)
      } else {
        return(
          tibble::tibble_row(
            estimator = estimator, 
            order = 2,
            entropy = the_entropy
          )
        )  
      }     
    } else {
      ## Tsallis entropy if estimator != "Lande" ----
      return(
        ent_tsallis(
          abd,
          q = 2, 
          estimator = estimator,
          level = NULL, 
          probability_estimator = probability_estimator,
          unveiling = unveiling,
          richness_estimator = richness_estimator,
          jack_alpha  = jack_alpha, 
          jack_max = jack_max,
          as_numeric = as_numeric,
          check_arguments = FALSE
        )
      )
    }
  }

  # Entropy at a level ----
  # If level is coverage, get size
  if (level < 1) {
    level <- coverage_to_size.numeric(
      abd, 
      sample_coverage = level,
      as_numeric = TRUE,
      check_arguments = FALSE
    )
  }
  if (level == sample_size) {
    # No interpolation/extrapolation needed: return the observed entropy
    the_entropy <- 1 - sum((abd / sample_size)^2)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "Sample", 
          order = 2,
          level = level,
          entropy = the_entropy
        )
      )  
    }
  }
  ## Valid interpolation and extrapolation ----
  the_entropy <- 1 - 1 / level - 
    (1 - 1 / level) * sum(abd * (abd - 1)) / sample_size / (sample_size - 1)
  if (as_numeric) {
    return(the_entropy)
  } else {
    return(
      tibble::tibble_row(
        estimator = "Chao2014", 
        order = 2,
        level = level,
        entropy = the_entropy
      )
    )  
  }
}


#' @rdname ent_simpson
#'
#' @export
ent_simpson.species_distribution <- function(
    x,
    estimator = c("Lande",
                  "UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive",
                  "Bonachela", "Holste"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    gamma = FALSE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  if (gamma) {
    return(
      ent_gamma.species_distribution(
        distribution = x,
        q = 2,
        estimator = estimator,
        level = level, 
        probability_estimator = probability_estimator,
        unveiling = unveiling,
        richness_estimator = richness_estimator,
        jack_alpha  = jack_alpha,
        jack_max = jack_max,
        as_numeric = FALSE
      )
    )
  } else {
    # Apply ent_simpson.numeric() to each site
    ent_simpson_list <- apply(
      # Eliminate site and weight columns
      x[, !(colnames(x) %in% non_species_columns)], 
      # Apply to each row
      MARGIN = 1,
      FUN = ent_simpson.numeric,
      # Arguments
      estimator = estimator,
      level = level, 
      probability_estimator = probability_estimator,
      unveiling = unveiling,
      richness_estimator = richness_estimator,
      jack_alpha  = jack_alpha, 
      jack_max = jack_max, 
      as_numeric = FALSE,
      check_arguments = FALSE
    )
    return(
      # Make a tibble with site, estimator and richness
      tibble::tibble(
        # Restore non-species columns
        x[colnames(x) %in% non_species_columns],
        # Coerce the list returned by apply into a dataframe
        do.call(rbind.data.frame, ent_simpson_list)
      )
    )
  }
}
