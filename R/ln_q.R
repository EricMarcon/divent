#' Deformed logarithm
#' 
#' Calculate the deformed logarithm of order *q*.
#' 
#' The deformed logarithm \insertCite{Tsallis1994}{divent} is defined as 
#' \eqn{\ln_q{x}=\frac{(x^{(1-q)}-1)}{(1-q)}}.
#' 
#' The shape of the deformed logarithm is similar to that of the regular one.
#' \eqn{\ln_1{x}=\log{x}}.
#'  
#' For \eqn{q>1}, \eqn{\ln_q{(+\infty)}=\frac{1}{(q-1)}}.
#'
#' @param x A numeric vector or array.
#' @param q A number.
#'
#' @return A vector of the same length as `x` containing the transformed values.
#' @export
#'  
#' @references
#' \insertAllCited{}
#'
#' @examples
#' curve(ln_q( 1/ x, q = 0), 0, 1, lty = 2, ylab = "Logarithm", ylim = c(0, 10))
#' curve(log(1 / x), 0, 1, lty = 1, n =1E4, add = TRUE)
#' curve(ln_q(1 / x, q = 2), 0, 1, lty = 3, add = TRUE)
#' legend("topright", 
#'   legend = c(
#'     expression(ln[0](1/x)),
#'     expression(log(1/x)),
#'     expression(ln[2](1/x))
#'   ),
#'   lty = c(2, 1, 3), 
#'   inset = 0.02
#'  )
#' 
ln_q <- function(x, q) {
  # Explicit recycling
  if (length(x) > length(q)) {
    q <- rep_len(q, length(x))
  } else if (length(x) < length(q)) {
    x <- rep_len(x, length(q))
  }
  
  # General formula
  the_ln_q <- (x^(1 - q) - 1) / (1 - q)
  # returns NaN if q==1: calculate log(x)
  the_ln_q[q == 1] <- log(x)[q == 1]
  # Force NaN for negative x
  the_ln_q[x < 0] <- NaN
  
  return (the_ln_q)
}
