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
#' @inheritParams check_divent_args
#' @param x An object of class [species_distribution] that contains several
#' communities or a matrix of abundances with communities in rows and species
#' in columns.
#' @param name The name of the metacommunity
#' @param ... Unused.
#'
#' @returns An object of class [abundances] or [probabilities] with a single
#' row, or a named vector if `as_numeric = TRUE`.
#'
#' @examples
#' metacommunity(paracou_6_abd)
#'
#' @name metacommunity
NULL


#' @rdname metacommunity
#'
#' @export
metacommunity <- function(
    x,
    name = "metacommunity",
    ...) {
  UseMethod("metacommunity")
}


#' @rdname metacommunity
#' @export
#'
#'
metacommunity.matrix <- function(
    x,
    name = "metacommunity",
    weights = rep(1, nrow(x)),
    as_numeric = TRUE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) {
    check_divent_args()
    if (length(weights) != nrow(x)) {
      cli::cli_abort("The length of 'weights' must be the number of communities")
    }
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }

  # Sample size
  sample_size <- sum(x)
  # Multiply abundances by weights and normalize so that
  # sample_size is the sum of sample sizes
  abd <- weights %*% x * sample_size / as.numeric(weights %*% rowSums(x))

  # Species names
  if (is.null(colnames(abd))) {
    colnames(abd) <- paste(
      "sp",
      formatC(
        seq_along(weights),
        width = ceiling(log10(length(weights))),
        flag = "0"
      ),
      sep = "_"
    )
  }

  if (as_numeric) {
    # Return a named vector
    return(abd[1, ])
  } else {
    # Build the tibble
    the_metacommunity <- tibble::as_tibble(
      cbind(
        data.frame(site = name, weight = sum(weights)),
        as.data.frame(abd)
      )
    )
    # Classes
    class(the_metacommunity) <- c("abundances", "species_distribution", class(the_metacommunity))
    return(the_metacommunity)
  }
}


#' @rdname metacommunity
#' @export
#'
metacommunity.species_distribution <- function(
    x,
    name = "metacommunity",
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  if (check_arguments) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species abundances must be positive.")
    }
  }

  # Select species columns
  species_columns <- !colnames(x) %in% non_species_columns
  # Extract abundances
  species_abd_prob <- as.matrix(x[, species_columns])
  # Sample size
  sample_size <- sum(species_abd_prob)
  # Aggregate abundances:
  # Multiply abundances by weights and normalize so that
  # sample_size is the sum of sample sizes
  abd_prob <- x$weight %*% species_abd_prob * sample_size /
    as.numeric(x$weight %*% rowSums(species_abd_prob))
  # Probabilities must still be divided by the number of communities
  if (is_probabilities(x)) {
    abd_prob <- abd_prob / nrow(species_abd_prob)
  }

  if (as_numeric) {
    # Return a named vector (abd_prob is a 1-row matrix)
    return(abd_prob[1, ])
  } else {
    # Build the tibble
    the_metacommunity <- tibble::as_tibble(
      cbind(
        data.frame(site = name, weight = sum(x$weight)),
        as.data.frame(abd)
      )
    )
    # Classes
    class(the_metacommunity) <- class(x)
    return(the_metacommunity)
  }
}

