# Shannon's Entropy of a Community

Estimate the entropy (Shannon 1948) of species from abundance or
probability data. Several estimators are available to deal with
incomplete sampling.

## Usage

``` r
ent_shannon(x, ...)

# S3 method for class 'numeric'
ent_shannon(
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

# S3 method for class 'species_distribution'
ent_shannon(
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

Bias correction requires the number of individuals.

See
[div_hill](https://ericmarcon.github.io/divent/dev/reference/div_hill.md)
for non-specific estimators. Shannon-specific estimators are from Miller
(1955) , Grassberger (2003) , Schürmann (2004) and Zhang (2012) . More
estimators can be found in the **entropy** package.

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
  
Grassberger P (2003). “Entropy Estimates from Insufficient Samplings.”
*arXiv Physics e-prints*, **0307138**(v2).  
  
Miller GA (1955). “Note on the Bias of Information Estimates.” In
Quastler H (ed.), *Information Theory in Psychology: Problems and
Methods*, 95–100. Free Press, Glencoe, Ill.  
  
Schürmann T (2004). “Bias Analysis in Entropy Estimation.” *Journal of
Physics A: Mathematical and General*, **37**(27), L295–L301.
[doi:10.1088/0305-4470/37/27/L02](https://doi.org/10.1088/0305-4470/37/27/L02)
.  
  
Shannon CE (1948). “A Mathematical Theory of Communication.” *The Bell
System Technical Journal*, **27**(3), 379–423, 623–656.
[doi:10.1002/j.1538-7305.1948.tb01338.x](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x)
.  
  
Zhang Z (2012). “Entropy Estimation in Turing's Perspective.” *Neural
Computation*, **24**(5), 1368–1389.
[doi:10.1162/NECO_a_00266](https://doi.org/10.1162/NECO_a_00266) .

## Examples

``` r
# Entropy of each community
ent_shannon(paracou_6_abd)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       1    4.57
#> 2 subplot_2   1.56 UnveilJ       1    4.73
#> 3 subplot_3   1.56 UnveilJ       1    4.65
#> 4 subplot_4   1.56 UnveilJ       1    4.55
# gamma entropy
ent_shannon(paracou_6_abd, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity UnveilJ       1    4.71

# At 80% coverage
ent_shannon(paracou_6_abd, level = 0.8)
#> # A tibble: 4 × 6
#>   site      weight estimator     order level entropy
#>   <chr>      <dbl> <chr>         <dbl> <dbl>   <dbl>
#> 1 subplot_1   1.56 Interpolation     1   304    4.10
#> 2 subplot_2   1.56 Interpolation     1   347    4.27
#> 3 subplot_3   1.56 Interpolation     1   333    4.23
#> 4 subplot_4   1.56 Interpolation     1   303    4.10
```
