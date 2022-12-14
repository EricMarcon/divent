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
#' @examples
#' curve(exp_q(x, q = 0), from = -5, to = 0, lty = 2)
#' curve(exp(x), from = -5, to = 0, lty= 1, add = TRUE)
#' curve(exp_q(x, q = 2), from = -5, to = 0, lty = 3, add = TRUE)
#' legend("topleft", 
#'   legend = c("exp_0(x)", "exp(x)", "exp_2(x)"), 
#'   lty = c(2, 1, 3), 
#'   inset = 0.02
#' )
#'
#' @references
#' \insertAllCited{}
#' 
exp_q <- function(x, q) {
  if (q == 1) {
    return (exp(x))
  } else {
    the_exp_q <- (x * (1 - q) + 1)^(1 / (1 - q))
    if (q > 1) {
      the_exp_q[x > 1 / (q - 1)] <- NA
    }
    return(the_exp_q)
  }
}
