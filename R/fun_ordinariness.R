#' Functional ordinariness of a community
#' 
#' The ordinariness of a species is the average similarity of its individuals 
#' with others \insertCite{Leinster2012}{divent}.
#' 
#' All species of the `species_distribution` must be found in the matrix of 
#' `similarities` if it is named.
#' If it is not, its size must equal the number of species.
#' Then, the order of species is assumed to be the same as that of the
#' `species_distribution`.
#'
#' @inheritParams check_divent_args
#'
#' @return A tibble with the ordinariness of each species, or a matrix if
#' argument `as_numeric` is `TRUE`.
#' @export
#'
#' @references
#' \insertAllCited{}
#' 
#' @examples
#' fun_ordinariness(paracou_6_abd, fun_similarity(paracou_6_fundist))
#' 
#' # Compare with probabilities
#' probabilities(paracou_6_abd)
#' # Decrease similarities so that ordinariness is close to probability
#' fun_ordinariness(paracou_6_abd, fun_similarity(paracou_6_fundist, rate = 100))
#' 
fun_ordinariness <- function (
    species_distribution,
    similarities = diag(
      sum(!colnames(species_distribution) %in% non_species_columns)
    ),
    as_numeric = FALSE,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) {
    check_divent_args()
    # Check species names
    similarities <- checked_matrix(similarities, species_distribution)
  }

  # Calculate species probabilities
  is_species_column <- !colnames(species_distribution) %in% non_species_columns
  prob <- probabilities(species_distribution)[, is_species_column]
  
  # Calculate ordinariness
  the_ordinariness <- as.matrix(prob) %*% t(similarities)
  
  if (as_numeric) {
    return(the_ordinariness)
  } else {
    return(
      tibble::tibble(
        species_distribution[, !is_species_column],
        as.data.frame(the_ordinariness)
      )
    )
  }
}
