# Diversity partition

Calculate \\\gamma\\, \\\beta\\ and \\\alpha\\ diversities of a
metacommunity.

## Usage

``` r
div_part(
  abundances,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Holste",
    "Marcon", "UnveilC", "UnveiliC", "ZhangGrabchak"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q_threshold = 10,
  check_arguments = TRUE
)
```

## Arguments

- abundances:

  an object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md).

- q:

  a number: the order of diversity.

- estimator:

  An estimator of diversity.

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

- q_threshold:

  the value of `q` above which diversity is computed directly with the
  naive estimator \\(\sum{p_s^q}^{\frac{1}{(1-q)}}\\, without computing
  entropy. When `q` is great, the exponential of entropy goes to
  \\0^{\frac{1}{(1-q)}}\\, causing rounding errors while the naive
  estimator of diversity is less and less biased.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A tibble with diversity values at each scale.

## Details

The function computes \\\gamma\\ diversity after building a
metacommunity from local communities according to their weight (Marcon
et al. 2014) . \\\alpha\\ entropy is the weighted mean local entropy,
converted into Hill numbers to obtain \\\alpha\\ diversity. \\\beta\\
diversity is obtained as the ratio of \\\gamma\\ to \\\alpha\\.

## References

Marcon E, Scotti I, Hérault B, Rossi V, Lang G (2014). “Generalization
of the Partitioning of Shannon Diversity.” *Plos One*, **9**(3), e90289.
[doi:10.1371/journal.pone.0090289](https://doi.org/10.1371/journal.pone.0090289)
.

## Examples

``` r
div_part(paracou_6_abd)
#> # A tibble: 7 × 6
#>   site          scale     estimator order diversity weight
#>   <chr>         <chr>     <chr>     <dbl>     <dbl>  <dbl>
#> 1 Metacommunity gamma     "UnveilJ"     1    111.     6.25
#> 2 Metacommunity beta      ""            1      1.09  NA   
#> 3 Metacommunity alpha     ""            1    102.    NA   
#> 4 subplot_1     community "UnveilJ"     1     96.3    1.56
#> 5 subplot_2     community "UnveilJ"     1    113.     1.56
#> 6 subplot_3     community "UnveilJ"     1    105.     1.56
#> 7 subplot_4     community "UnveilJ"     1     94.6    1.56
```
