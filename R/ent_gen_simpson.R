#' Generalized Simpson's Entropy
#'
#' Estimate the Generalized Simpson's entropy of species from abundance
#' or probability data.
#'
#' The Generalized Simpson's Entropy \insertCite{Zhang2010}{divent} of order \eqn{k} is,
#' in the species accumulation curve,the probability for the individual sampled
#' in rank \eqn{k + 1} to belong to a new species.
#' It is a measure of diversity so long as \eqn{k} is lower than the number
#' of species \insertCite{Grabchak2016}{divent}.
#'
#' Bias correction requires the number of individuals.
#' It is limited to orders \eqn{r} less than or equal to the number of individuals
#' in the community \insertCite{Zhang2014}{divent}.
#'
#' Generalized Simpson's diversity cannot be estimated at a specified level
#' of interpolation or extrapolation, and diversity partitioning is not available.
#'
#' @note The unbiased estimator is calculated by the [EntropyEstimation::GenSimp.z]
#' function of the **EntropyEstimation** package.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated entropy.
#'
#' @seealso [div_gen_simpson]
#'
#' #' @references
#' \insertAllCited{}
#'
#' @examples
#' # Entropy of each community
#' ent_gen_simpson(paracou_6_abd, k = 50)
#' # gamma entropy
#' ent_gen_simpson(paracou_6_abd, k = 50, gamma = TRUE)
#'
#' @name ent_gen_simpson
NULL


#' @rdname ent_gen_simpson
#'
#' @export
ent_gen_simpson <- function(x, ...) {
  UseMethod("ent_gen_simpson")
}


#' @rdname ent_gen_simpson
#'
#' @param estimator An estimator of entropy.
#'
#' @export
ent_gen_simpson.numeric <- function(
    x,
    k = 1,
    estimator = c("Zhang", "naive"),
    as_numeric  = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }

  # Entropy of a vector of probabilities ----
  if (abs(sum(x) - 1) < length(x) * .Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    the_entropy <- sum(x * (1 - x)^k)
    if (as_numeric) {
      return(the_entropy)
    } else {
      return(
        tibble::tibble_row(
          estimator = "naive",
          order = k,
          entropy = the_entropy
        )
      )
    }
  }

  # Integer values needed
  if (!is_integer_values(x)) {
    cli::cli_abort("Species abundances must be integers.")
  }
  # Computed by the EntropyEstimation package
  the_entropy <- EntropyEstimation::GenSimp.z(x = x, r = k)
  # Return the entropy
  if (as_numeric) {
    return(the_entropy)
  } else {
    return(
      tibble::tibble_row(
        estimator = estimator,
        order = k,
        entropy = the_entropy
      )
    )
  }
}


#' @rdname ent_gen_simpson
#'
#' @export
ent_gen_simpson.species_distribution <- function(
    x,
    k = 1,
    estimator = c("Zhang", "naive"),
    gamma = FALSE,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }

  if (gamma) {
    # Build the metacommunity
    abd <- metacommunity.abundances(x, as_numeric = TRUE, check_arguments = FALSE)
    if (estimator != "naive" & !is_integer_values(abd)) {
      cli::cli_abort(
        paste(
          "The weights of communities yield non-integer abundances in",
          "the metacommunity that only allow using the naive estimator."
        )
      )
    }
    return(
      ent_gen_simpson.numeric(
        x = abd,
        k = k,
        estimator = estimator,
        as_numeric = as_numeric,
        check_arguments = FALSE
      )
    )
  } else {
    # Apply ent_gen_simpson.numeric() to each site
    ent_simpson_sites <- apply(
      # Eliminate site and weight columns
      x[, !colnames(x) %in% non_species_columns],
      # Apply to each row
      MARGIN = 1,
      FUN = ent_gen_simpson.numeric,
      # Arguments
      k = k,
      estimator = estimator,
      as_numeric = as_numeric,
      check_arguments = FALSE
    )
    if (as_numeric) {
      return(ent_simpson_sites)
    } else {
      return(
        # Make a tibble with site, estimator and richness
        tibble::tibble(
          # Restore non-species columns
          x[colnames(x) %in% non_species_columns],
          # Coerce the list returned by apply into a dataframe
          do.call(rbind.data.frame, ent_simpson_sites)
        )
      )
    }
  }
}
