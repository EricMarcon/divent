#' Log-Series Distribution
#' 
#' Random generation for the log-series distribution.
#' 
#' Fast implementation of the random generation of a log-series distribution
#' \insertCite{Fisher1943}{divent}.
#' 
#' The complete set of functions (including density, distribution function and quantiles)
#' can be found in package _sads_ but this implementation of the random generation is much faster.
#' 
#' If `size` is too large, i.e. `size` + 1 can't be distinguished from `size` due to rounding,
#' then an error is raised.
#'
#' @inheritParams check_divent_args
#' @param n Number of observations.
#' @param size The size of the distribution.
#' @param alpha Fisher's \eqn{\alpha}.
#'
#' @returns A numeric vector with the random values drawn from the log-series distribution.
#' @export
#'
#' @references
#' \insertAllCited{}

#' @examples
#' # Generate a community made of 10000 individuals with alpha=40
#' size <- 1E4
#' alpha <- 40
#' species_number <- -alpha * log(alpha / (size + alpha))
#' abundances <- rlseries(species_number, size = 1E5, alpha = 40)
#' # rcommunity() may be a better choice here
#' autoplot(rcommunity(1, size = 1E4, alpha = 40, distribution = "lseries"))
rlseries <- function(
    n, 
    size, 
    alpha, 
    show_progress = TRUE, 
    check_arguments = TRUE) {
  # adapted from Dan Lunn, http://www.stats.ox.ac.uk/~dlunn/BS1_05/BS1_Rcode.pdf
  if (size + 1 == size || size - 1 == size) {
    stop("size is too large to simulate the distribution.")
  }
  # Fisher's x
  x <- size / (size + alpha)
  # Draw random numbers between 0 and 1 that are quantiles of the cumulative function
  u <- sort(stats::runif(n))
  # Prepare the corresponding number of abundances
  abd <- numeric(n)
  # Prepare a progress bar
  if (show_progress & interactive()) {
    cli::cli_progress_bar(
      total = n,
      format = "{cli::pb_spin} Drawing log-series values {cli::pb_current}/{cli::pb_total}"
    )
  }
  # Number of drawn values +1
  next_value <- 1
  # k is the abundance
  k <- 1
  # Calculate the probability at k=1
  P <- -x / log(1 - x)
  # Store it in the cumulative function
  F <- P
  # Repeat while all values are not drawn
  while (next_value <= n) {
    if (F > u[next_value]) {
      # Retain k as the next value
      abd[next_value] <- k
      # Increment the next value
      next_value <- next_value + 1
      # Update the progress bar
      if (show_progress & interactive()) cli::cli_progress_update()
    } else {
      # Probability at k+1 obtained from that at k
      P <- P * k * x / (k + 1)
      # Increment the cumulated probability
      F <- F + P
      # Increment k
      k <- k + 1
    }
  }
  return(abd)
}
