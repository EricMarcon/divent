#' Aggregate communities into a metacommunity
#'
#' @param abundances An object of class [abundances] that contains several communities.
#' @param name The name of the metacommunity
#' @param as_numeric If `TRUE`, a number is returned rather than a tibble.
#' @param check_arguments If `TRUE`, the function arguments are verified.
#' Should be set to `FALSE` to save time when the arguments have been checked elsewhere.
#'
#' @return An object of class [abundances] with a single row.
#' @export
#'
#' @examples
#' metacommunity(paracou_6_abd)
metacommunity <- function(
    abundances, 
    name = "metacommunity",
    as_numeric = FALSE,
    check_arguments = TRUE) {
  
  if (check_arguments) check_divent_args()
  
  # Select species columns
  species_columns <- !(colnames(abundances) %in% c("site", "weight"))
  # Extract abundances
  species_abd <- as.matrix(abundances[, species_columns])
  sample_size <- sum(species_abd)
  # Multiply them by weights and normalize so that 
  # sample_size is the sum of sample sizes
  abd <- abundances$weight %*% species_abd *
    sample_size / as.numeric(abundances$weight %*% rowSums(species_abd))
  if (as_numeric) {
    return(as.numeric(abd))
  } else {
    # Build the tibble
    the_metacommunity <- tibble::as_tibble(
      cbind(
        data.frame(site = name, weight = sample_size),
        as.data.frame(abd)
      )
    )
    # Restore exact species names (spaces may have been transformed into "_")
    colnames(the_metacommunity[, species_columns]) <- colnames(abundances[, species_columns])
    # Classes
    class(the_metacommunity) <- c("abundances", "species_distribution", class(the_metacommunity))
    return(the_metacommunity)
  }
}
