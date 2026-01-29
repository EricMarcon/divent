# Similarity-Based Diversity Profile of a Community

Calculate the diversity profile of a community, i.e. its
similarity-based diversity against its order.

## Usage

``` r
profile_similarity(
  x,
  similarities,
  orders = seq(from = 0, to = 2, by = 0.1),
  ...
)

# S3 method for class 'numeric'
profile_similarity(
  x,
  similarities = diag(length(x)),
  orders = seq(from = 0, to = 2, by = 0.1),
  estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang", "UnveilC", "UnveiliC",
    "naive"),
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  sample_coverage = NULL,
  as_numeric = FALSE,
  n_simulations = 0,
  alpha = 0.05,
  bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
profile_similarity(
  x,
  similarities = diag(sum(!colnames(x) %in% non_species_columns)),
  orders = seq(from = 0, to = 2, by = 0.1),
  estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang", "UnveilC", "UnveiliC",
    "naive"),
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  gamma = FALSE,
  n_simulations = 0,
  alpha = 0.05,
  bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object, that may be a numeric vector containing abundances or
  probabilities, or an object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md).

- similarities:

  a similarity matrix, that can be obtained by
  [fun_similarity](https://ericmarcon.github.io/divent/reference/fun_similarity.md).
  Its default value is the identity matrix.

- orders:

  The orders of diversity used to build the profile.

- ...:

  Unused.

- estimator:

  An estimator of entropy.

- probability_estimator:

  a string containing one of the possible estimators of the probability
  distribution (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only for extrapolation.

- unveiling:

  a string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only for extrapolation.

- richness_estimator:

  an estimator of richness to evaluate the total number of species, see
  [div_richness](https://ericmarcon.github.io/divent/reference/div_richness.md).
  used for interpolation and extrapolation.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/reference/coverage.md).

- sample_coverage:

  the sample coverage of `x` calculated elsewhere. Used to calculate the
  gamma diversity of meta-communities, see details.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- n_simulations:

  The number of simulations used to estimate the confidence envelope of
  the profile.

- alpha:

  The risk level, 5% by default, of the confidence envelope of the
  profile.

- bootstrap:

  the method used to obtain the probabilities to generate bootstrapped
  communities from observed abundances. If "Marcon2012", the
  probabilities are simply the abundances divided by the total number of
  individuals (Marcon et al. 2012) . If "Chao2013" or "Chao2015" (by
  default), a more sophisticated approach is used (see
  [as_probabilities](https://ericmarcon.github.io/divent/reference/species_distribution.md))
  following Chao et al. (2013) or Chao and Jost (2015) .

- show_progress:

  if TRUE, a progress bar is shown during long computations.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- gamma:

  if `TRUE`, \\\gamma\\ diversity, i.e. diversity of the metacommunity,
  is computed.

## Value

A tibble with the site names, the estimators used and the estimated
diversity at each order. This is an object of class "profile" that can
be plotted.

## Details

A bootstrap confidence interval can be produced by simulating
communities (their number is `n_simulations`) with
[rcommunity](https://ericmarcon.github.io/divent/reference/rcommunity.md)
and calculating their profiles. Simulating communities implies a
downward bias in the estimation: rare species of the actual community
may have abundance zero in simulated communities. Simulated diversity
values are recentered so that their mean is that of the actual
community.

## References

Chao A, Jost L (2015). “Estimating Diversity and Entropy Profiles via
Discovery Rates of New Species.” *Methods in Ecology and Evolution*,
**6**(8), 873–882.
[doi:10.1111/2041-210X.12349](https://doi.org/10.1111/2041-210X.12349)
.  
  
Chao A, Wang Y, Jost L (2013). “Entropy and the Species Accumulation
Curve: A Novel Entropy Estimator via Discovery Rates of New Species.”
*Methods in Ecology and Evolution*, **4**(11), 1091–1100.
[doi:10.1111/2041-210x.12108](https://doi.org/10.1111/2041-210x.12108)
.  
  
Marcon E, Hérault B, Baraloto C, Lang G (2012). “The Decomposition of
Shannon's Entropy and a Confidence Interval for *Beta* Diversity.”
*Oikos*, **121**(4), 516–522.
[doi:10.1111/j.1600-0706.2011.19267.x](https://doi.org/10.1111/j.1600-0706.2011.19267.x)
.

## Examples

``` r
# Similarity matrix
Z <- fun_similarity(paracou_6_fundist)
# Profile
profile_similarity(paracou_6_abd, similarities = Z, q = 2)
#> # A tibble: 84 × 4
#>    site      estimator order diversity
#>    <chr>     <chr>     <dbl>     <dbl>
#>  1 subplot_1 UnveilJ     0        1.31
#>  2 subplot_1 UnveilJ     0.1      1.31
#>  3 subplot_1 UnveilJ     0.2      1.31
#>  4 subplot_1 UnveilJ     0.3      1.31
#>  5 subplot_1 UnveilJ     0.4      1.31
#>  6 subplot_1 UnveilJ     0.5      1.31
#>  7 subplot_1 UnveilJ     0.6      1.31
#>  8 subplot_1 UnveilJ     0.7      1.31
#>  9 subplot_1 UnveilJ     0.8      1.31
#> 10 subplot_1 UnveilJ     0.9      1.31
#> # ℹ 74 more rows
```
