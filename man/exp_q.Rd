% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exp_q.R
\name{exp_q}
\alias{exp_q}
\title{Deformed exponential}
\usage{
exp_q(x, q)
}
\arguments{
\item{x}{A numeric vector or array.}

\item{q}{A number.}
}
\value{
A vector of the same length as \code{x} containing the transformed values.
}
\description{
Calculate the deformed exponential of order \emph{q}.
}
\details{
The deformed exponential is the reciprocal of the deformed logarithm
\insertCite{Tsallis1994}{divent}, see \link{ln_q}.
It is defined as \eqn{(x(1-q)+1)^{\frac{1}{(1-q)}}}.
For \eqn{q>1}, \eqn{\ln_q{(+\infty)}=\frac{1}{(q-1)}}
so \eqn{\exp_q{(x)}} is not defined for \eqn{x>\frac{1}{(q-1)}}.
}
\examples{
curve(exp_q(x, q = 0), from = -5, to = 0, lty = 2)
curve(exp(x), from = -5, to = 0, lty= 1, add = TRUE)
curve(exp_q(x, q = 2), from = -5, to = 0, lty = 3, add = TRUE)
legend("bottomright", 
  legend = c(
    expression(exp[0](x)),
    expression(exp(x)),
    expression(exp[2](x))
  ),
  lty = c(2, 1, 3), 
  inset = 0.02
)

}
\references{
\insertAllCited{}
}
