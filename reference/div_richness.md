# Number of Species of a Community

Estimate the number of species from abundance or probability data.
Several estimators are available to deal with incomplete sampling.

## Usage

``` r
div_richness(x, ...)

# S3 method for class 'numeric'
div_richness(
  x,
  estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
div_richness(
  x,
  estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
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

  Unused. The metacommunity if built by combining the community
  abundances with respect to their weight.

- estimator:

  An estimator of richness to evaluate the total number of species.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- level:

  The level of interpolation or extrapolation. It may be a sample size
  (an integer) or a sample coverage (a number between 0 and 1). The
  asymptotic `estimator` is used in extrapolation (i.e. a `level`
  greater than the sample size).

- probability_estimator:

  A string containing one of the possible estimators of the probability
  distribution (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only by the estimator of richness "rarefy".

- unveiling:

  A string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only by the estimator of richness "rarefy".

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
numbers of species.

## Details

Bias correction requires the number of individuals. Chao's estimation
techniques are from Chao et al. (2014) and Chiu et al. (2014) . The
Jackknife estimator is calculated by a straight adaptation of the code
by Ji-Ping Wang (jackknife in package **SPECIES**). The optimal order is
selected according to Burnham and Overton (1978); Burnham and Overton
(1979) . Many other estimators are available elsewhere, the ones
implemented here are necessary for other entropy estimations.

Richness can be estimated at a specified `level` of interpolation or
extrapolation, either a chosen sample size or sample coverage (Chiu et
al. 2014) , rather than its asymptotic value. Extrapolation relies on
the estimation of the asymptotic richness. If `probability_estimator` is
"naive", then the asymptotic estimation of richness is made using the
chosen `estimator`, else the asymptotic distribution of the community is
derived and its estimated richness adjusted so that the richness of a
sample of this distribution of the size of the actual sample has the
richness of the actual sample.

## References

Burnham KP, Overton WS (1978). “Estimation of the Size of a Closed
Population When Capture Probabilities Vary among Animals.” *Biometrika*,
**65**(3), 625–633.
[doi:10.2307/2335915](https://doi.org/10.2307/2335915) .  
  
Burnham KP, Overton WS (1979). “Robust Estimation of Population Size
When Capture Probabilities Vary among Animals.” *Ecology*, **60**(5),
927–936. [doi:10.2307/1936861](https://doi.org/10.2307/1936861) .  
  
Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .  
  
Chiu C, Wang Y, Walther BA, Chao A (2014). “An Improved Nonparametric
Lower Bound of Species Richness via a Modified Good-Turing Frequency
Formula.” *Biometrics*, **70**(3), 671–682.
[doi:10.1111/biom.12200](https://doi.org/10.1111/biom.12200) , 24945937.

## Examples

``` r
# Diversity of each community
div_richness(paracou_6_abd)
#> # A tibble: 4 × 5
#>   site      weight estimator   order diversity
#>   <chr>      <dbl> <chr>       <dbl>     <dbl>
#> 1 subplot_1   1.56 Jackknife 3     0       355
#> 2 subplot_2   1.56 Jackknife 2     0       348
#> 3 subplot_3   1.56 Jackknife 2     0       315
#> 4 subplot_4   1.56 Jackknife 2     0       296
# gamma diversity
div_richness(paracou_6_abd, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator   order diversity
#>   <chr>         <chr>       <dbl>     <dbl>
#> 1 Metacommunity Jackknife 2     0       473

# At 80% coverage
div_richness(paracou_6_abd, level = 0.8)
#> # A tibble: 4 × 6
#>   site      weight estimator order level diversity
#>   <chr>      <dbl> <chr>     <dbl> <dbl>     <dbl>
#> 1 subplot_1   1.56 SAC           0   304      106.
#> 2 subplot_2   1.56 SAC           0   347      125.
#> 3 subplot_3   1.56 SAC           0   333      117.
#> 4 subplot_4   1.56 SAC           0   303      109.
```
