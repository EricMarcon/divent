# Generalized Simpson's Diversity

Estimate the diversity sensu stricto, i.e. the effective number of
species (Grabchak et al. 2017) from abundance or probability data.

## Usage

``` r
div_gen_simpson(x, k = 1, ...)

# S3 method for class 'numeric'
div_gen_simpson(
  x,
  k = 1,
  estimator = c("Zhang", "naive"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
div_gen_simpson(
  x,
  k = 1,
  estimator = c("Zhang", "naive"),
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

- k:

  the order of Hurlbert's diversity.

- ...:

  Unused.

- estimator:

  An estimator of asymptotic diversity.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A tibble with the site names, the estimators used and the estimated
diversity.

## Details

Bias correction requires the number of individuals.

Estimation techniques are from Zhang and Grabchak (2016) . It is limited
to orders \\k\\ less than or equal to the number of individuals in the
community.

Generalized Simpson's diversity cannot be estimated at a specified level
of interpolation or extrapolation, and diversity partitioning is not
available.

## References

Grabchak M, Marcon E, Lang G, Zhang Z (2017). “The Generalized Simpson's
Entropy Is a Measure of Biodiversity.” *Plos One*, **12**(3), e0173305.
[doi:10.1371/journal.pone.0173305](https://doi.org/10.1371/journal.pone.0173305)
.  
  
Zhang Z, Grabchak M (2016). “Entropic Representation and Estimation of
Diversity Indices.” *Journal of Nonparametric Statistics*, **28**(3),
563–575.
[doi:10.1080/10485252.2016.1190357](https://doi.org/10.1080/10485252.2016.1190357)
.

## See also

[ent_gen_simpson](https://ericmarcon.github.io/divent/dev/reference/ent_gen_simpson.md)

## Examples

``` r
# Diversity of each community
div_gen_simpson(paracou_6_abd, k = 50)
#> # A tibble: 4 × 5
#>   site      weight estimator order diversity
#>   <chr>      <dbl> <chr>     <dbl>     <dbl>
#> 1 subplot_1   1.56 Zhang        50      1.01
#> 2 subplot_2   1.56 Zhang        50      1.02
#> 3 subplot_3   1.56 Zhang        50      1.01
#> 4 subplot_4   1.56 Zhang        50      1.01
```
