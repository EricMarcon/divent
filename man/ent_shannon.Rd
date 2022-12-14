% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ent_shannon.R
\name{ent_shannon}
\alias{ent_shannon}
\alias{ent_shannon.numeric}
\alias{ent_shannon.species_distribution}
\title{Shannon's Entropy of a Community}
\usage{
ent_shannon(x, ...)

\method{ent_shannon}{numeric}(
  x,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Grassberger2003",
    "Holste", "Miller", "Schurmann", "ZhangHz"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

\method{ent_shannon}{species_distribution}(
  x,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Grassberger2003",
    "Holste", "Miller", "Schurmann", "ZhangHz"),
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
or an object of class \link{abundances} or \link{probabilities}.}

\item{...}{Unused.}

\item{estimator}{An estimator of entropy.}

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

\item{as_numeric}{If \code{TRUE}, a number or a numeric vector is returned rather than a tibble.}

\item{check_arguments}{If \code{TRUE}, the function arguments are verified.
Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.}

\item{gamma}{If \code{TRUE}, \eqn{\gamma} diversity, i.e. diversity of the metacommunity, is computed.}
}
\value{
A tibble with the site names, the estimators used and the estimated entropy.
}
\description{
Estimate the entropy \insertCite{Shannon1948}{divent} of species from abundance
or probability data.
Several estimators are available to deal with incomplete sampling.
}
\details{
Bias correction requires the number of individuals.

See \link{div_hill} for non-specific estimators.
Shannon-specific estimators are from \insertCite{Miller1955;textual}{divent},
\insertCite{Grassberger2003;textual}{divent}, \insertCite{Schurmann2004;textual}{divent}
and \insertCite{Zhang2012;textual}{divent}.
More estimators can be found in the \strong{entropy} package.

Entropy can be estimated at a specified level of interpolation or
extrapolation, either a chosen sample size or sample coverage
\insertCite{Chao2014}{divent}, rather than its asymptotic value.
See \link{accum_tsallis} for details.
}
\examples{
# Entropy of each community
ent_shannon(paracou_6_abd)
# gamma entropy
ent_shannon(paracou_6_abd, gamma = TRUE)

# At 80\% coverage
ent_shannon(paracou_6_abd, level = 0.8)

}
\references{
\insertAllCited{}
}
