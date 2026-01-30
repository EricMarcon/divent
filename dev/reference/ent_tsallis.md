# Tsallis Entropy of a Community

Estimate the entropy of species from abundance or probability data.
Several estimators are available to deal with incomplete sampling.

## Usage

``` r
ent_tsallis(x, q = 1, ...)

# S3 method for class 'numeric'
ent_tsallis(
  x,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  sample_coverage = NULL,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
ent_tsallis(
  x,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  gamma = FALSE,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object, that may be a numeric vector containing abundances or
  probabilities, or an object of class
  [abundances](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md).

- q:

  a number: the order of diversity.

- ...:

  Unused.

- estimator:

  An estimator of entropy.

- level:

  the level of interpolation or extrapolation. It may be a sample size
  (an integer) or a sample coverage (a number between 0 and 1). If not
  `NULL`, the asymptotic `estimator` is ignored.

- probability_estimator:

  a string containing one of the possible estimators of the probability
  distribution (see
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md)).
  Used only for extrapolation.

- unveiling:

  a string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species (see
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md)).
  Used only for extrapolation.

- richness_estimator:

  an estimator of richness to evaluate the total number of species, see
  [div_richness](https://ericmarcon.github.io/divent/dev/reference/div_richness.md).
  used for interpolation and extrapolation.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/dev/reference/coverage.md).

- sample_coverage:

  the sample coverage of `x` calculated elsewhere. Used to calculate the
  gamma diversity of meta-communities, see details.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- gamma:

  if `TRUE`, \\\gamma\\ diversity, i.e. diversity of the metacommunity,
  is computed.

## Value

A tibble with the site names, the estimators used and the estimated
entropy.

## Details

Bias correction requires the number of individuals. See
[div_hill](https://ericmarcon.github.io/divent/dev/reference/div_hill.md)
for estimators.

Entropy can be estimated at a specified level of interpolation or
extrapolation, either a chosen sample size or sample coverage (Chao et
al. 2014) , rather than its asymptotic value. See
[accum_tsallis](https://ericmarcon.github.io/divent/dev/reference/accum_hill.md)
for details.

## References

Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .

## Examples

``` r
# Entropy of each community
ent_tsallis(paracou_6_abd, q = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       2   0.976
#> 2 subplot_2   1.56 UnveilJ       2   0.978
#> 3 subplot_3   1.56 UnveilJ       2   0.980
#> 4 subplot_4   1.56 UnveilJ       2   0.972
# gamma entropy
ent_tsallis(paracou_6_abd, q = 2, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity UnveilJ       2   0.979

# At 80% coverage
ent_tsallis(paracou_6_abd, level = 0.8)
#> # A tibble: 4 × 6
#>   site      weight estimator     order level entropy
#>   <chr>      <dbl> <chr>         <dbl> <dbl>   <dbl>
#> 1 subplot_1   1.56 Interpolation     1   304    4.10
#> 2 subplot_2   1.56 Interpolation     1   347    4.27
#> 3 subplot_3   1.56 Interpolation     1   333    4.23
#> 4 subplot_4   1.56 Interpolation     1   303    4.10
```
