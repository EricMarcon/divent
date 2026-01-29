# Abundance Frequency Count of a Community

Count the number of species observed the same number of times.

## Usage

``` r
abd_freq_count(
  abd,
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  check_arguments = TRUE
)
```

## Arguments

- abd:

  A numeric vector containing species abundances.

- level:

  the level of interpolation or extrapolation. If `NULL` (by default),
  observed data is simply used and all other arguments are ignored.

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

  A string containing an estimator recognized by
  [div_richness](https://ericmarcon.github.io/divent/reference/div_richness.md)
  to evaluate the total number of species in
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md).
  Used only for extrapolation.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/reference/coverage.md).

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A two-column tibble. The first column contains the number of
observations, the second one the number of species observed this number
of times.

## Details

The Abundance Frequency Count (Chao and Jost 2015) is the number of
species observed each number of times. It is a way to summarize the
species distribution.

It can be estimated at a specified level of interpolation or
extrapolation. Extrapolation relies on the estimation of the estimation
of the asymptotic distribution of the community by
[probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)
and eq. (5) of (Chao et al. 2014) .

## References

Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .  
  
Chao A, Jost L (2015). “Estimating Diversity and Entropy Profiles via
Discovery Rates of New Species.” *Methods in Ecology and Evolution*,
**6**(8), 873–882.
[doi:10.1111/2041-210X.12349](https://doi.org/10.1111/2041-210X.12349) .

## Examples

``` r
abd_freq_count(as.numeric(paracou_6_abd[1, ]))
#> # A tibble: 25 × 2
#>    abundance number_of_species
#>        <dbl>         <int[1d]>
#>  1         1                84
#>  2         2                35
#>  3         3                19
#>  4         4                14
#>  5         5                 5
#>  6         6                 4
#>  7         7                 3
#>  8         8                 3
#>  9        10                 2
#> 10        11                 1
#> # ℹ 15 more rows
```
