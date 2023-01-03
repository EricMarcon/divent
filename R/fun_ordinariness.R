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
  similarities <- similarities_checked(similarities, species_distribution)
  
  # Calculate species probabilities
  is_species_column <- !(colnames(species_distribution) %in% non_species_columns)
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


#' Check Similarity matrix
#' 
#' Verify that a similarity matrix fits a species distribution and filter it
#' so that its elements are the same, in the same order.
#'
#' @param similarities A similarity matrix
#' @param species_distribution A species distribution, or a named vector.
#'
#' @return A similarity matrix that corresponds to the species distribution.
#' @noRd
#'
similarities_checked <- function(
    similarities,
    species_distribution) {
  
  # No names needed
  if (is.null(colnames(similarities))) {
    # Similarities may not be named
    if (ncol(similarities) != length(species_names)) {
      stop("If the similarity matrix is not named, then its size must fit the number of species.")
    } else {
      # Do not change the similarities
      return(similarities)
    }
  }
  
  # Get species names
  if (is_species_distribution(species_distribution)) {
    is_species_column <- !(colnames(species_distribution) %in% non_species_columns)
    species_names <- colnames(species_distribution)[is_species_column]
  } else if (is.vector(species_distribution)) {
    species_names <- names(species_distribution)
  }

  # Stop if some species are not in the matrix of similarities
  if (length(species_names) == 0) stop("There are no species in the distribution")
  if (length(setdiff(species_names, colnames(similarities))) != 0) {
    stop("Some species are missing in the similarity matrix.")    
  } 
  
  # Filter and reorder the similarity matrix
  return(similarities[species_names, species_names])
}
