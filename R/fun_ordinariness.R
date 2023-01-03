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
      sum(!(colnames(species_distribution) %in% non_species_columns))
    ),
    as_numeric = FALSE,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  
  # Check species names
  is_species_column <- !(colnames(species_distribution) %in% non_species_columns)
  species_names <- colnames(species_distribution)[is_species_column]
  # Similarities may not be named
  if (is.null(colnames(similarities))) {
    if (ncol(similarities) != length(species_names)) {
      stop("If the similarity matrix is not named, then its size must fit the number of species.")
    }
  } else {
    if (length(setdiff(species_names, colnames(similarities))) != 0) {
      stop("Some species are missing in the similarity matrix.")    
    } else {
      # Filter and reorder the similarity matrix
      similarities <- similarities[species_names, species_names]
    }
  }
  
  # Calculate species probabilities
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
