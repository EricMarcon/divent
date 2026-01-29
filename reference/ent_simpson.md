# Simpson's Entropy of a Community

Estimate the entropy (Simpson 1949) of species from abundance or
probability data. Several estimators are available to deal with
incomplete sampling.

## Usage

``` r
ent_simpson(x, ...)

# S3 method for class 'numeric'
ent_simpson(
  x,
  estimator = c("Lande", "UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger",
    "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
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

# S3 method for class 'species_distribution'
ent_simpson(
  x,
  estimator = c("Lande", "UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger",
    "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
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
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md).

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
[div_hill](https://ericmarcon.github.io/divent/reference/div_hill.md)
for non-specific estimators.

Simpson-specific estimator is from Lande (1996) .

Entropy can be estimated at a specified level of interpolation or
extrapolation, either a chosen sample size or sample coverage (Chao et
al. 2014) , rather than its asymptotic value. See
[accum_tsallis](https://ericmarcon.github.io/divent/reference/accum_hill.md)
for details.

## References

Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .  
  
Lande R (1996). “Statistics and Partitioning of Species Diversity, and
Similarity among Multiple Communities.” *Oikos*, **76**(1), 5–13.
[doi:10.2307/3545743](https://doi.org/10.2307/3545743) .  
  
Simpson EH (1949). “Measurement of Diversity.” *Nature*, **163**(4148),
688. [doi:10.1038/163688a0](https://doi.org/10.1038/163688a0) .

## Examples

``` r
# Entropy of each community
ent_simpson(paracou_6_abd)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 Lande         2   0.976
#> 2 subplot_2   1.56 Lande         2   0.978
#> 3 subplot_3   1.56 Lande         2   0.980
#> 4 subplot_4   1.56 Lande         2   0.972
# gamma entropy
ent_simpson(paracou_6_abd, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity Lande         2   0.979

# At 80% coverage
ent_simpson(paracou_6_abd, level = 0.8)
#> # A tibble: 4 × 6
#>   site      weight estimator order level entropy
#>   <chr>      <dbl> <chr>     <dbl> <dbl>   <dbl>
#> 1 subplot_1   1.56 Chao2014      2   304   0.973
#> 2 subplot_2   1.56 Chao2014      2   347   0.975
#> 3 subplot_3   1.56 Chao2014      2   333   0.977
#> 4 subplot_4   1.56 Chao2014      2   303   0.969
```
