#' Functional ordinariness of a community
#' 
#' The ordinariness of a species is the average similarity of its individuals 
#' with others \insertCite{Leinster2012}{divent}.
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
fun_ordinariness <- function (
    species_distribution,
    similarities,
    as_numeric = FALSE,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  
  # Species names
  is_species_column <- !(colnames(species_distribution) %in% non_species_columns)
  species_names <- colnames(species_distribution)[is_species_column]
  if (length(setdiff(species_names, colnames(similarities))) != 0) {
    stop("Some species are missing in the similarity matrix.")    
  } else {
    # Filter and reorder the similarity matrix
    similarities <- similarities[species_names, species_names]
  }
  
  # Calculate species probabilities
  prob <- probabilities(species_distribution)[, is_species_column]
  
  # Calculate ordinariness
  the_similarities <- as.matrix(prob) %*% t(similarities)
  
  if (as_numeric) {
    return(the_similarities)
  } else {
    return(
      tibble::tibble(
        species_distribution[, !is_species_column],
        as.data.frame(the_similarities)
      )
    )
  }
}
