# Sample Coverage of a Community

`coverage()` calculates an estimator of the sample coverage of a
community described by its abundance vector. `coverage_to_size()`
estimates the sample size corresponding to the chosen sample coverage.

## Usage

``` r
coverage(x, ...)

# S3 method for class 'numeric'
coverage(
  x,
  estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  level = NULL,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'abundances'
coverage(
  x,
  estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  level = NULL,
  ...,
  check_arguments = TRUE
)

coverage_to_size(x, ...)

# S3 method for class 'numeric'
coverage_to_size(
  x,
  sample_coverage,
  estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'abundances'
coverage_to_size(
  x,
  sample_coverage,
  estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object.

- ...:

  Unused.

- estimator:

  An estimator of the sample coverage. "ZhangHuang" is the most accurate
  but does not allow choosing a `level`. "Good"'s estimator only allows
  interpolation, i.e. estimation of the coverage of a subsample. "Chao"
  allows estimation at any `level`, including extrapolation. "Turing" is
  the simplest estimator, equal to 1 minus the number of singletons
  divided by the sample size.

- level:

  The level of interpolation or extrapolation, i.e. an abundance.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- sample_coverage:

  The target sample coverage.

## Value

`coverage()` returns a named number equal to the calculated sample
coverage. The name is that of the estimator used.

`coverage_to_size()` returns a number equal to the sample size
corresponding to the chosen sample coverage.

## Details

The sample coverage \\C\\ of a community is the total probability of
occurrence of the species observed in the sample. \\1-C\\ is the
probability for an individual of the whole community to belong to a
species that has not been sampled.

The historical estimator is due to Turing (Good 1953) . It only relies
on singletons (species observed only once). Chao's (Chao and Shen 2010)
estimator uses doubletons too and Zhang-Huang's (Chao et al. 1988; Zhang
and Huang 2007) uses the whole distribution.

If `level` is not `NULL`, the sample coverage is interpolated or
extrapolated. Interpolation by the Good estimator relies on the equality
between sampling deficit and the generalized Simpson entropy (Good 1953)
. The Chao et al. (2014) estimator allows extrapolation, reliable up a
level equal to the double size of the sample.

## References

Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .  
  
Chao A, Lee S, Chen T (1988). “A Generalized Good's Nonparametric
Coverage Estimator.” *Chinese Journal of Mathematics*, **16**, 189–199.
43836340.  
  
Chao A, Shen T (2010). *Program SPADE: Species Prediction and Diversity
Estimation. Program and User's Guide.*. CARE.  
  
Good IJ (1953). “The Population Frequency of Species and the Estimation
of Population Parameters.” *Biometrika*, **40**(3/4), 237–264.
[doi:10.1093/biomet/40.3-4.237](https://doi.org/10.1093/biomet/40.3-4.237)
.  
  
Zhang Z, Huang H (2007). “Turing's Formula Revisited.” *Journal of
Quantitative Linguistics*, **14**(2-3), 222–241.
[doi:10.1080/09296170701514189](https://doi.org/10.1080/09296170701514189)
.

## Examples

``` r
coverage(paracou_6_abd)
#> # A tibble: 4 × 4
#>   site      weight estimator  coverage
#>   <chr>      <dbl> <chr>         <dbl>
#> 1 subplot_1   1.56 ZhangHuang    0.911
#> 2 subplot_2   1.56 ZhangHuang    0.893
#> 3 subplot_3   1.56 ZhangHuang    0.912
#> 4 subplot_4   1.56 ZhangHuang    0.902
coverage_to_size(paracou_6_abd, sample_coverage = 0.9)
#> # A tibble: 4 × 4
#>   site      weight sample_coverage  size
#>   <chr>      <dbl>           <dbl> <dbl>
#> 1 subplot_1   1.56             0.9   819
#> 2 subplot_2   1.56             0.9   940
#> 3 subplot_3   1.56             0.9   826
#> 4 subplot_4   1.56             0.9   778
```
