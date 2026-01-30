# Hurlbert Entropy of a Community

Estimate the Hurlbert entropy (Hurlbert 1971) of species from abundance
or probability data. Several estimators are available to deal with
incomplete sampling.

## Usage

``` r
ent_hurlbert(x, k = 2, ...)

# S3 method for class 'numeric'
ent_hurlbert(
  x,
  k = 2,
  estimator = c("Hurlbert", "naive"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
ent_hurlbert(
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
  [abundances](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md).

- k:

  the order of Hurlbert's diversity.

- ...:

  Unused.

- estimator:

  An estimator of entropy.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A tibble with the site names, the estimators used and the estimated
entropy.

## Details

Bias correction requires the number of individuals. See
[div_hurlbert](https://ericmarcon.github.io/divent/dev/reference/div_hurlbert.md)
for estimators.

Hurlbert's entropy cannot be estimated at a specified level of
interpolation or extrapolation, and entropy partitioning is not
available.

## References

Hurlbert SH (1971). “The Nonconcept of Species Diversity: A Critique and
Alternative Parameters.” *Ecology*, **52**(4), 577–586.
[doi:10.2307/1934145](https://doi.org/10.2307/1934145) .

## See also

[div_hurlbert](https://ericmarcon.github.io/divent/dev/reference/div_hurlbert.md)

## Examples

``` r
# Entropy of each community
ent_hurlbert(paracou_6_abd, k = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 Hurlbert      2    1.98
#> 2 subplot_2   1.56 Hurlbert      2    1.98
#> 3 subplot_3   1.56 Hurlbert      2    1.98
#> 4 subplot_4   1.56 Hurlbert      2    1.97
```
