% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlseries.R
\name{rlseries}
\alias{rlseries}
\title{Log-Series Distribution}
\usage{
rlseries(n, size, alpha, show_progress = TRUE, check_arguments = TRUE)
}
\arguments{
\item{n}{Number of observations.}

\item{size}{The size of the distribution.}

\item{alpha}{Fisher's \eqn{\alpha}.}

\item{show_progress}{If TRUE, a progress bar is shown during long computations.}

\item{check_arguments}{If \code{TRUE}, the function arguments are verified.
Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.}
}
\value{
A numeric vector with the random values drawn from the log-series distribution.
}
\description{
Random generation for the log-series distribution.
}
\details{
Fast implementation of the random generation of a log-series distribution
\insertCite{Fisher1943}{divent}.

The complete set of functions (including density, distribution function and quantiles)
can be found in package \emph{sads} but this implementation of the random generation is much faster.

If \code{size} is too large, i.e. \code{size} + 1 can't be distinguished from \code{size} due to rounding,
then an error is raised.
}
\examples{
# Generate a community made of 10000 individuals with alpha=40
size <- 1E4
alpha <- 40
species_number <- -alpha * log(alpha / (size + alpha))
abundances <- rlseries(species_number, size = 1E5, alpha = 40)
# rCommunity() may be a better choice
autoplot(rcommunity(n = 1, size = 1E4, alpha = 40, distribution = "lseries"))
}
\references{
\insertAllCited{}
}
