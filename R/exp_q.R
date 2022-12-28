#' Deformed exponential
#' 
#' Calculate the deformed exponential of order *q*.
#' 
#' The deformed exponential is the reciprocal of the deformed logarithm
#' \insertCite{Tsallis1994}{divent}, see [ln_q].
#' It is defined as \eqn{(x(1-q)+1)^{\frac{1}{(1-q)}}}.
#' For \eqn{q>1}, \eqn{\ln_q{(+\infty)}=\frac{1}{(q-1)}} 
#' so \eqn{\exp_q{(x)}} is not defined for \eqn{x>\frac{1}{(q-1)}}.
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
#' curve(exp_q(x, q = 0), from = -5, to = 0, lty = 2)
#' curve(exp(x), from = -5, to = 0, lty= 1, add = TRUE)
#' curve(exp_q(x, q = 2), from = -5, to = 0, lty = 3, add = TRUE)
#' legend("bottomright", 
#'   legend = c(
#'     expression(exp[0](x)),
#'     expression(exp(x)),
#'     expression(exp[2](x))
#'   ),
#'   lty = c(2, 1, 3), 
#'   inset = 0.02
#' )
#' 
exp_q <- function(x, q) {
  # Explicit recycling
  if (length(x) > length(q)) {
    q <- rep_len(q, length(x))
  } else if (length(x) < length(q)) {
    x <- rep_len(x, length(q))
  }
  
  # General formula
  the_exp_q <- (x * (1 - q) + 1)^(1 / (1 - q))
  # returns 1 if q==1: calculate exp(x)
  the_exp_q[q == 1] <- exp(x)[q == 1]
  # Force NaN for x out of limits
  the_exp_q[q > 1 & x > 1 / (q - 1)] <- NaN
  
  return(the_exp_q)
}
