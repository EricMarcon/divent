#' Grassberger's expectation of n^q
#' 
#' Expected value of \eqn{n^q} when \eqn{n} follows a Poisson distribution 
#' of parameter \eqn{n}.
#' 
#' The expectation of \eqn{n^q} when \eqn{n} follows a Poisson distribution 
#' was derived by \insertCite{Grassberger1988;textual}{divent}.
#' 
#' The function is computed using the [beta] function.
#' Its value is 0 for \eqn{n-q+1<0}.
#'
#' @param n A positive integer vector.
#' @param q A positive number.
#'
#' @return A vector of the same length as n containing the transformed values.
#' @export
#'
#' @examples
#' # Compare
#' n <- c(2, 3)
#' e_n_q(n, q = 2)
#' # with
#' n^2
#' 
#' # Result is 1
#' e_n_q(n, q = 0)
#' # Result is 0
#' e_n_q(n, q=5)
#' @references
#' \insertAllCited{}
e_n_q <- function(n, q) {
  if (q == 0) {
    return (rep(1, length(n)))
  } else {
    # beta cannot be computed for n - q + 1 < 0 (so warnings must be suppressed) 
    # but the value is 0 then
    beta_value  <- suppressWarnings(gamma(q) / beta(n - q + 1, q))
    beta_value[n - q + 1 < 0] <- 0
    # (-1)^n is problematic for long vectors (returns NA for large values). 
    # It is replaced by 1 - n %% 2 * 2 (n is rounded if is not an integer)
    the_e_n_q <- beta_value - 
      (1 - round(n) %% 2 * 2) * gamma(1 + q) * sin(pi * q) / pi / (n + 1)   
    return (the_e_n_q)
  }
}                                
