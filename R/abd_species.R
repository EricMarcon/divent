#' Abundances of Communities
#' 
#' Utilities for community abundances (objects of class "abundances").
#' 
#' `abd_species()` returns a tibble containing the species abundance columns only,
#' to simplify numeric operations.
#' `prob_species()` returns the same tibble but values are probabilities.
#' 
#' `abd_sum()` returns the sample sizes of the communities.
#'
#' @inheritParams check_divent_args
#'
#' @examples
#' abd_species(paracou_6_abd)
#' abd_sum(paracou_6_abd)
#' 
#' @name abd_species
NULL


#' @rdname abd_species
#'
#' @export
abd_species <- function(
    abundances, 
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  
  return(
    abundances[, !colnames(abundances) %in% non_species_columns]
  )
}


#' @rdname abd_species
#'
#' @export
abd_sum <- function(
    abundances, 
    as_numeric = FALSE, 
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  
  the_abd_sum <- rowSums(
    abundances[, !colnames(abundances) %in% non_species_columns]
  )
  if (as_numeric) {
    return(the_abd_sum)
  } else {
    return(
      # Make a tibble with site etc., estimator and richness
      tibble::tibble(
        # Restore non-species columns
        abundances[colnames(abundances) %in% non_species_columns],
        abundance = the_abd_sum
      )
    )
  }
}


#' @rdname abd_species
#'
#' @export
prob_species <- function(
    species_distribution, 
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  
  abundances <- species_distribution[
    , 
    !colnames(species_distribution) %in% non_species_columns
  ]
  sample_sizes <- rowSums(abundances)
  # Divide each column by sample_sizes
  return(
    dplyr::mutate(
      abundances,
      dplyr::across(
        dplyr::everything(),
        ~ .x / sample_sizes
      )
    )
  )
}
