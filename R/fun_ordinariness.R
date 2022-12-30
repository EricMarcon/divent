#' Functional ordinariness of a community
#' 
#' TODO
#'
#' @inheritParams check_divent_args
#' @param distribution TODO
#' @param similarities TODO
#'
#' @return A tibble with the ordinariness of each species.
#' @export
#'
#' @examples
#' fun_ordinariness(paracou_6_abd, fun_similarity(paracou_6_fundist))
#' 
fun_ordinariness <- function (
    distribution,
    similarities,
    as_numeric = FALSE,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  
  # Species names
  is_species_column <- !(colnames(distribution) %in% non_species_columns)
  species_names <- colnames(distribution)[is_species_column]
  if (length(setdiff(species_names, colnames(similarities))) != 0) {
    stop("Some species are missing in the similarity matrix.")    
  } else {
    # Filter and reorder the similarity matrix
    similarities <- similarities[species_names, species_names]
  }
  
  # Calculate species probabilities
  prob <- probabilities(distribution)[, is_species_column]
  
  # Calculate ordinariness
  the_similarities <- as.matrix(prob) %*% t(similarities)
  
  if (as_numeric) {
    return(the_similarities)
  } else {
    return(
      tibble::tibble(
        distribution[, !is_species_column],
        as.data.frame(the_similarities)
      )
    )
  }
}
