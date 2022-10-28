#' Tsallis Entropy of a Community
#' 
#' Estimate the entropy species from abundance or probability data.
#' Several estimators are available to deal with incomplete sampling.
#' 
#' Bias correction requires the number of individuals. 
#'
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances]  or [probabilities]
#' @param ... Unused.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#'
#' @return A tibble with the site names, the estimators used and the estimated entropy.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' # TODO
#' # ent_tsallis(paracou_6_abd, q = 2)
#' 
#' @name ent_tsallis
NULL


#' @rdname ent_tsallis
#'
#' @export
ent_tsallis <- function(x, q = 1, ...) {
  UseMethod("ent_tsallis")
}


#' @rdname ent_tsallis
#'
#' @param q The order or diversity.
#' @param estimator An estimator of entropy. 
#' @param level The level of interpolation or extrapolation. 
#' It may be a chosen sample size (an integer) or a sample coverage 
#' (a number between 0 and 1). 
#' Richness extrapolation require its asymptotic estimation depending on the 
#' choice of the `estimator`.
#' @param probability_estimator One of the estimators of a probability distribution: 
#' "naive" (the default value), or "Chao2013", "Chao2015", "ChaoShen" to estimate
#' the probabilities of the observed species in the asymptotic distribution.
#' @param unveiling One of the possible unveiling methods to estimate the probabilities 
#' of the unobserved species: "none" (default, no species is added), "uniform" 
#' (all unobserved species have the same probability) or "geometric" (the 
#' unobserved species distribution is geometric).
#' @param richness_estimator An estimator of richness to evaluate the total number of species,
#' see [div_richness].
#' @param jack_alpha The risk level, 5% by default, used to optimize the jackknife order.
#' @param jack_max The highest jackknife order allowed. Default is 10. 
#' 
#' @export
ent_tsallis.numeric <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak"),
    level = NULL, 
    probability_estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    sample_coverage = NULL,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    prob <- x[x > 0]
    ent_species <- prob * ln_q(1 / prob, q = q)
    return (entropy)
    return(
      tibble::tibble_row(
        estimator = "naive", 
        richness = sum(ent_species)
      )
    )
  }
  
  # Eliminate 0
  abd <- x[x > 0]
  # Sample size
  sample_size <- sum(abd)
  # Number of observed species
  s_obs <- length(abd)
  
  # Entropy of a vector of abundances ----
  if (is.null(level)) {
    ## Exit if x contains no or a single species ----
    if (length(abd) < 2) {
      if (length(abd) == 0) {
        return(tibble::tibble_row(estimator = "No Species", entropy = NA))
      } else {
        return(tibble::tibble_row(estimator = "Single Species", entropy = 1))
      }
    } else {
      # Probabilities instead of abundances
      if (sample_size < 2) {
        warning("Entropy estimators can't apply to probability data. Estimator forced to 'naive'")
        estimator <- "naive"
      }
    }
    
    ## Metacommunity estimation ----
    # abd may not be integers, sample_coverage is given
    # sample_coverage is between 0 and 1 (by check_arguments), sum(abd) must be an integer.
    # estimator may be ChaoShen or Marcon (max(ChoShen, Grassberger))
    if (
        !is.null(sample_coverage) & 
        is_integer_values(sample_size) & 
        (estimator == "ChaoShen" | estimator == "Marcon")) {
      
      cp <- sample_coverage * abd / sample_size
      chao_shen <- sum(cp^q * ln_q(1 / cp, q = q) /(1 - (1 - cp)^sample_size))
      if (estimator == "Marcon") {
        # Calculate Grassberger's correction
        if (q == 1) {
          grassberger <- sum(
            abd / sample_size * (log(sample_size) - digamma(abd) - 
            (1 - round(abd) %% 2 * 2) / (abd + 1))
          )
        } else {
          grassberger <- (1 - sample_size^(-q) * sum(e_n_q(abd, q = q))) / (q - 1)
        }
      } else grassberger <- 0
      # Take the max
      if (chao_shen > grassberger) {
        return(tibble::tibble_row(estimator = "ChaoShen", entropy = chao_shen))
      } else {
        return(tibble::tibble_row(estimator = "Grassberger", entropy = grassberger))
      }
    }
  }  
  
}


#' @rdname ent_tsallis
#'
#' @export
ent_tsallis.abundances <- function(
    x, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak"),
    level = NULL, 
    probability_estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
    unveiling = c("none", "uniform", "geometric"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(x < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Apply ent_tsallis.numeric() to each site
  ent_tsallis_list <- apply(
    # Eliminate site and weight columns
    x[, !(colnames(x) %in% c("site", "weight"))], 
    # Apply to each row
    MARGIN = 1,
    FUN = ent_tsallis.numeric,
    # Arguments
    q = q,
    estimator = estimator,
    level = level, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max, 
    check_arguments = FALSE
  )
  
  return(
    # Make a tibble with site, estimator and richness
    tibble::tibble(
      # Do not assume column site exists
      x[colnames(x) == "site"],
      # Coerce the list returned by apply into a dataframe
      do.call(rbind.data.frame, ent_tsallis_list)
    )
  )
  
}
