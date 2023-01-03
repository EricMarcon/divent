#' Aggregate communities into a metacommunity
#' 
#' Abundances of communities are summed according to their weights to obtain
#' the abundances of the metacommunity.
#' 
#' The total abundance of the metacommunity is by design equal to the sum 
#' of community abundances so that the information used by diversity estimators.
#' A consequence is that equal weights lead to a metacommunity whose species 
#' abundances are the sum of community species abundances.
#' 
#' If community weights are not equal then the metacommunity abundances are 
#' in general not integer.
#' Most diversity estimators can't be applied to non-integer abundances but the knowledge
#' of the sample coverage of each community allow "ChaoShen" and "Grassberger"
#' estimators.
#' 
#'
#' @inheritParams check_divent_args
#' @param abundances An object of class [abundances] that contains several communities.
#' @param name The name of the metacommunity
#'
#' @return An object of class [abundances] with a single row.
#' @export
#'
#' @examples
#' metacommunity(paracou_6_abd)
#' 
metacommunity <- function(
    abundances, 
    name = "metacommunity",
    as_numeric = FALSE,
    check_arguments = TRUE) {
  
  if (any(check_arguments)) check_divent_args()
  
  # Select species columns
  species_columns <- !colnames(abundances) %in% non_species_columns
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
