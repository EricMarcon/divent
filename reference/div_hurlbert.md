# Hurlbert Diversity of a Community

Estimate the diversity sensu stricto, i.e. the effective number of
species Dauby and Hardy (2012) from abundance or probability data.

## Usage

``` r
div_hurlbert(x, k = 1, ...)

# S3 method for class 'numeric'
div_hurlbert(
  x,
  k = 2,
  estimator = c("Hurlbert", "naive"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
div_hurlbert(
  x,
  k = 2,
  estimator = c("Hurlbert", "naive"),
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

Several estimators are available to deal with incomplete sampling.

Bias correction requires the number of individuals.

Estimation techniques are from Hurlbert (1971) .

Hurlbert's diversity cannot be estimated at a specified level of
interpolation or extrapolation, and diversity partioning is not
available.

## References

Dauby G, Hardy OJ (2012). “Sampled-Based Estimation of Diversity Sensu
Stricto by Transforming Hurlbert Diversities into Effective Number of
Species.” *Ecography*, **35**(7), 661–672.
[doi:10.1111/j.1600-0587.2011.06860.x](https://doi.org/10.1111/j.1600-0587.2011.06860.x)
.  
  
Hurlbert SH (1971). “The Nonconcept of Species Diversity: A Critique and
Alternative Parameters.” *Ecology*, **52**(4), 577–586.
[doi:10.2307/1934145](https://doi.org/10.2307/1934145) .

## See also

[ent_hurlbert](https://ericmarcon.github.io/divent/reference/ent_hurlbert.md)

## Examples

``` r
# Diversity of each community
div_hurlbert(paracou_6_abd, k = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator order diversity
#>   <chr>      <dbl> <chr>     <dbl>     <dbl>
#> 1 subplot_1   1.56 Hurlbert      2      42.3
#> 2 subplot_2   1.56 Hurlbert      2      44.6
#> 3 subplot_3   1.56 Hurlbert      2      48.9
#> 4 subplot_4   1.56 Hurlbert      2      36.0
```
