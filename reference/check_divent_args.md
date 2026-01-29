# check_divent_args

Checks the arguments of a function of the package divent

## Usage

``` r
check_divent_args(
  abd = NULL,
  abundances = NULL,
  alpha = NULL,
  as_numeric = NULL,
  bootstrap = NULL,
  check_arguments = NULL,
  correction = NULL,
  coverage_estimator = NULL,
  distances = NULL,
  distribution = NULL,
  estimator = NULL,
  fisher_alpha = NULL,
  gamma = NULL,
  global = NULL,
  jack_alpha = NULL,
  jack_max = NULL,
  k = NULL,
  level = NULL,
  n = NULL,
  n_simulations = NULL,
  normalize = NULL,
  orders = NULL,
  prob = NULL,
  prob_geom = NULL,
  probability_estimator = NULL,
  q = NULL,
  q_threshold = NULL,
  r = NULL,
  rate = NULL,
  richness_estimator = NULL,
  sample_coverage = NULL,
  sd_lnorm = NULL,
  show_progress = NULL,
  similarities = NULL,
  size = NULL,
  species_number = NULL,
  species_distribution = NULL,
  thomas_mu = NULL,
  thomas_scale = NULL,
  tree = NULL,
  use.names = NULL,
  unveiling = NULL,
  weights = NULL,
  w_max = NULL,
  w_mean = NULL,
  w_min = NULL,
  weibull_scale = NULL,
  weibull_shape = NULL,
  X = NULL,
  win = NULL
)
```

## Arguments

- abd:

  a numeric vector containing abundances.

- abundances:

  an object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md).

- alpha:

  the risk level, 5% by default.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- bootstrap:

  the method used to obtain the probabilities to generate bootstrapped
  communities from observed abundances. If "Marcon2012", the
  probabilities are simply the abundances divided by the total number of
  individuals (Marcon et al. 2012) . If "Chao2013" or "Chao2015" (by
  default), a more sophisticated approach is used (see
  [as_probabilities](https://ericmarcon.github.io/divent/reference/species_distribution.md))
  following Chao et al. (2013) or Chao and Jost (2015) .

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- correction:

  the edge-effect correction to apply when estimating the number of
  neighbors.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/reference/coverage.md).

- distances:

  a distance matrix or an object of class
  [stats::dist](https://rdrr.io/r/stats/dist.html).

- distribution:

  The distribution of species abundances. May be "lnorm" (log-normal),
  "lseries" (log-series), "geom" (geometric) or "bstick" (broken stick).

- estimator:

  an estimator of asymptotic entropy, diversity or richness.

- fisher_alpha:

  Fisher's \\\alpha\\ in the log-series distribution.

- gamma:

  if `TRUE`, \\\gamma\\ diversity, i.e. diversity of the metacommunity,
  is computed.

- global:

  if `TRUE`, a global envelope sensu (Duranton and Overman 2005) is
  calculated.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- k:

  the order of Hurlbert's diversity.

- level:

  the level of interpolation or extrapolation. It may be a sample size
  (an integer) or a sample coverage (a number between 0 and 1). If not
  `NULL`, the asymptotic `estimator` is ignored.

- n:

  the number of observations.

- n_simulations:

  the number of simulations used to estimate the confidence envelope.

- normalize:

  if `TRUE`, phylogenetic is normalized: the height of the tree is set
  to 1.

- orders:

  The orders of diversity.

- prob:

  a numeric vector containing probabilities.

- prob_geom:

  the proportion of resources taken by successive species of the
  geometric distribution.

- probability_estimator:

  a string containing one of the possible estimators of the probability
  distribution (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only for extrapolation.

- q:

  a number: the order of diversity.

- q_threshold:

  the value of `q` above which diversity is computed directly with the
  naive estimator \\(\sum{p_s^q}^{\frac{1}{(1-q)}}\\, without computing
  entropy. When `q` is great, the exponential of entropy goes to
  \\0^{\frac{1}{(1-q)}}\\, causing rounding errors while the naive
  estimator of diversity is less and less biased.

- r:

  a vector of distances.

- rate:

  the decay rate of the exponential similarity.

- richness_estimator:

  an estimator of richness to evaluate the total number of species, see
  [div_richness](https://ericmarcon.github.io/divent/reference/div_richness.md).
  used for interpolation and extrapolation.

- sample_coverage:

  the sample coverage of `x` calculated elsewhere. Used to calculate the
  gamma diversity of meta-communities, see details.

- sd_lnorm:

  the simulated log-normal distribution standard deviation. This is the
  standard deviation on the log scale.

- show_progress:

  if TRUE, a progress bar is shown during long computations.

- similarities:

  a similarity matrix, that can be obtained by
  [fun_similarity](https://ericmarcon.github.io/divent/reference/fun_similarity.md).
  Its default value is the identity matrix.

- size:

  the number of individuals to draw in each community.

- species_number:

  the number of species.

- species_distribution:

  an object of class
  [species_distribution](https://ericmarcon.github.io/divent/reference/species_distribution.md).

- thomas_mu:

  in Thomas point patterns, the mean number of points per cluster. The
  intensity of the Poisson process of cluster centers is calculated as
  the number of points (`size`) per area divided by `thomas_mu`.

- thomas_scale:

  in Thomas point patterns, the standard deviation of random
  displacement (along each coordinate axis) of a point from its cluster
  center.

- tree:

  an ultrametric, phylogenetic tree. May be an object of class
  [phylo_divent](https://ericmarcon.github.io/divent/reference/phylo_divent.md),
  [ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html),
  [ade4::phylog](https://adeverse.github.io/ade4/reference/phylog.html)
  or [stats::hclust](https://rdrr.io/r/stats/hclust.html).

- use.names:

  if `TRUE`, the names of the `species_distribution` are kept in the
  matrix or vector they are converted to.

- unveiling:

  a string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only for extrapolation.

- weights:

  the weights of the sites of the species distributions.

- w_max:

  the maximum weight in a uniform distribution.

- w_mean:

  the mean weight in an exponential distribution (i.e. the negative of
  the inverse of the decay rate).

- w_min:

  the minimum weight in a uniform, exponential or Weibull distribution.

- weibull_scale:

  the scale parameter in a Weibull distribution.

- weibull_shape:

  the shape parameter in a Weibull distribution.

- X:

  a spatialized community (A
  [dbmss::wmppp](https://ericmarcon.github.io/dbmss/reference/wmppp.html)
  object with `PointType` values as species names.)

- win:

  the window containing the point pattern. It is an
  [spatstat.geom::owin](https://rdrr.io/pkg/spatstat.geom/man/owin.html)
  object. Default is a 1x1 square.

## Value

Returns `TRUE` or stops if a problem is detected.

## Details

The function compares the arguments passed to its parent function to the
type they should be and performs some extra tests, *e.g.* probabilities
must be positive and sum to 1. It stops if an argument is not correct.

The function is always called without arguments. Its arguments exist
only for documentation.
