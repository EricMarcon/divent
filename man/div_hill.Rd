% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/div_hill.R
\name{div_hill}
\alias{div_hill}
\alias{div_hill.numeric}
\alias{div_hill.species_distribution}
\title{Hill number of a Community}
\usage{
div_hill(x, q = 1, ...)

\method{div_hill}{numeric}(
  x,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  sample_coverage = NULL,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

\method{div_hill}{species_distribution}(
  x,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  gamma = FALSE,
  ...,
  check_arguments = TRUE
)
}
\arguments{
\item{x}{An object, that may be a numeric vector containing abundances or probabilities,
or an object of class \link{abundances}  or \link{probabilities}.}

\item{q}{The order of diversity.}

\item{...}{Unused.}

\item{estimator}{An estimator of asymptotic diversity.}

\item{level}{The level of interpolation or extrapolation.
It may be a sample size (an integer) or a sample coverage
(a number between 0 and 1).
If not \code{NULL}, the asymptotic \code{estimator} is ignored.}

\item{probability_estimator}{A string containing one of the possible estimators
of the probability distribution (see \link{probabilities}).
Used only for extrapolation.}

\item{unveiling}{A string containing one of the possible unveiling methods
to estimate the probabilities of the unobserved species (see \link{probabilities}).
Used only for extrapolation.}

\item{richness_estimator}{An estimator of richness to evaluate the total number of species,
see \link{div_richness}. Used for interpolation and extrapolation.}

\item{jack_alpha}{The risk level, 5\% by default, used to optimize the jackknife order.}

\item{jack_max}{The highest jackknife order allowed. Default is 10.}

\item{coverage_estimator}{An estimator of sample coverage used by \link{coverage}.}

\item{sample_coverage}{The sample coverage of \code{x} calculated elsewhere.
Used to calculate the gamma diversity of meta-communities, see details.}

\item{as_numeric}{If \code{TRUE}, a number or a numeric vector is returned rather than a tibble.}

\item{check_arguments}{If \code{TRUE}, the function arguments are verified.
Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.}

\item{gamma}{If \code{TRUE}, \eqn{\gamma} diversity, i.e. diversity of the metacommunity, is computed.}
}
\value{
A tibble with the site names, the estimators used and the estimated diversity.
}
\description{
Estimate the diversity sensu stricto, i.e. the \insertCite{Hill1973;textual}{divent}
number of species from abundance or probability data.
}
\details{
Several estimators are available to deal with incomplete sampling.

Bias correction requires the number of individuals.

Estimation techniques are from \insertCite{Chao2003;textual}{divent},
\insertCite{Grassberger1988;textual}{divent},\insertCite{Holste1998;textual}{divent},
\insertCite{Bonachela2008;textual}{divent}, \insertCite{Marcon2014a;textual}{divent}
which is actually the max value of "ChaoShen" and "Grassberger",
\insertCite{Zhang2014a;textual}{divent}, \insertCite{Chao2014c;textual}{divent},
\insertCite{Chao2015;textual}{divent} and \insertCite{Marcon2015a;textual}{divent}.

The \code{ChaoJost} estimator \insertCite{Chao2013,Chao2015;textual}{divent} contains
an unbiased part concerning observed species, equal to that of
\insertCite{Zhang2014a;textual}{divent}, and a (biased) estimator of the remaining
bias based on the estimation of the species-accumulation curve.
It is very efficient but slow if the number of individuals is more than a few hundreds.

The unveiled estimators rely on \insertCite{Chao2014c;textual}{divent},
completed by \insertCite{Marcon2015a;textual}{divent}.
The actual probabilities of observed species are estimated and completed by
a geometric distribution of the probabilities of unobserved species.
The number of unobserved species is estimated by the Chao1 estimator (\code{UnveilC}),
following \insertCite{Chao2014c;textual}{divent}, or by the iChao1 (\code{UnveiliC})
or the jackknife (\code{UnveilJ}).
The \code{UnveilJ} correction often has a lower bias but a greater variance
\insertCite{Marcon2015a}{divent}.
It is a good first choice thanks to the versatility of the jackknife
estimator of richness.

Estimators by \insertCite{Bonachela2008;textual}{divent} and
\insertCite{Holste1998;textual}{divent} are rarely used.

To estimate \eqn{\gamma} diversity, the size of a metacommunity (see
\link{metacommunity}) is unknown so it has to be set according to a rule which does
not ensure that its abundances are integer values.
Then, classical bias-correction methods do not apply.
Providing the \code{sample_coverage} argument allows applying the \code{ChaoShen} and
\code{Grassberger} corrections to estimate quite well the entropy.

Diversity can be estimated at a specified level of interpolation or
extrapolation, either a chosen sample size or sample coverage
\insertCite{Chao2014}{divent}, rather than its asymptotic value.
See \link{accum_hill} for details.
}
\examples{
# Diversity of each community
div_hill(paracou_6_abd, q = 2)
# gamma diversity
div_hill(paracou_6_abd, q = 2, gamma = TRUE)

# At 80\% coverage
div_hill(paracou_6_abd, q = 2, level = 0.8)

}
\references{
\insertAllCited{}
}
