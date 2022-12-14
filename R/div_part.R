#' Diversity partition
#' 
#' Calculate \eqn{\gamma}, \eqn{\beta} and \eqn{\alpha} diversities of a metacommunity.
#' 
#' The function computes \eqn{\gamma} diversity after building a metacommunity 
#' from local communities according to their weight \insertCite{Marcon2014a}{divent}.
#' \eqn{\alpha} entropy is the weighted mean local entropy, converted into Hill
#' numbers to obtain \eqn{\alpha} diversity.
#' \eqn{\beta} diversity is obtained as the ratio of \eqn{\gamma} to \eqn{\alpha}.
#'
#' @inheritParams check_divent_args
#' @param estimator An estimator of diversity.
#'
#' @return A tibble with diversity values at each scale.
#' @export
#' 
#' @references
#' \insertAllCited{}
#'
#' @examples
#' div_part(paracou_6_abd)
#' 
div_part <- function(
    abundances, 
    q = 1, 
    estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", 
                  "Holste", "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak"),
    level = NULL, 
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
    jack_alpha  = 0.05, 
    jack_max = 10,
    check_arguments = TRUE) {
 
  if (any(check_arguments)) check_divent_args()
  estimator <- match.arg(estimator) 
  probability_estimator <- match.arg(probability_estimator) 
  unveiling <- match.arg(unveiling) 
  richness_estimator <- match.arg(richness_estimator) 
  if (any(abundances < 0)) stop("Species probabilities or abundances must be positive.")
  
  # Gamma diversity ----
  # Calculate gamma entropy
  ent_gamma <- ent_gamma(
    x = abundances,
    q = q,
    estimator = estimator,
    level = level, 
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha,
    jack_max = jack_max
  )
  # Calculate diversity
  div_gamma <- dplyr::mutate(
    ent_gamma, 
    diversity = exp_q(.data$entropy, q = q),
    .keep = "unused"
  )
  # Add a scale column
  div_gamma <- dplyr::bind_cols(scale = "gamma", div_gamma)

  # Alpha diversity ----
  # Apply div_hill.numeric() to each site
  ent_sites <- apply(
    # Eliminate site and weight columns
    abundances[, !(colnames(abundances) %in% c("site", "weight"))], 
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
    as_numeric = TRUE,
    check_arguments = FALSE
  )
  # Calculate alpha entropy
  ent_alpha <- as.numeric(
    abundances$weight %*% ent_sites / sum(abundances$weight)
  )
  # Alpha diversity
  div_alpha <- tibble::tibble_row(
    scale = "alpha",
    estimator = "",
    order = q,
    diversity = exp_q(ent_alpha, q = q)
  )
  # Beta diversity
  div_beta <- tibble::tibble_row(
    scale = "beta",
    estimator = "",
    order = q,
    diversity = div_gamma$diversity / div_alpha$diversity
  )
  
  # Summarize
  div_metacommunity <- dplyr::bind_cols(
    site = "metacommunity",
    dplyr::bind_rows(div_gamma, div_beta, div_alpha)
  )
    
  # Site diversity
  div_sites <- div_hill.species_distribution(
    x = abundances,
    q = q,
    estimator = estimator,
    level = level,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    richness_estimator = richness_estimator,
    jack_alpha  = jack_alpha, 
    jack_max = jack_max, 
    gamma = FALSE,
    check_arguments = FALSE
  )
  
  return(
    dplyr::bind_rows(
      div_metacommunity, 
      dplyr::bind_cols(scale = "community", div_sites)
    )
  )
}
