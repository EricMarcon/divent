# Generalized Simpson's Entropy

Estimate the Generalized Simpson's entropy of species from abundance or
probability data.

## Usage

``` r
ent_gen_simpson(x, ...)

# S3 method for class 'numeric'
ent_gen_simpson(
  x,
  k = 1,
  estimator = c("Zhang", "naive"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
ent_gen_simpson(
  x,
  k = 1,
  estimator = c("Zhang", "naive"),
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

- k:

  the order of Hurlbert's diversity.

- estimator:

  An estimator of entropy.

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

The Generalized Simpson's Entropy (Zhang and Zhou 2010) of order \\k\\
is, in the species accumulation curve,the probability for the individual
sampled in rank \\k + 1\\ to belong to a new species. It is a measure of
diversity so long as \\k\\ is lower than the number of species (Grabchak
et al. 2017) .

Bias correction requires the number of individuals. It is limited to
orders \\r\\ less than or equal to the number of individuals in the
community (Zhang and Grabchak 2016) .

Generalized Simpson's diversity cannot be estimated at a specified level
of interpolation or extrapolation, and diversity partitioning is not
available.

## Note

The unbiased estimator is calculated by the
[EntropyEstimation::GenSimp.z](https://rdrr.io/pkg/EntropyEstimation/man/GenSimp.z.html)
function of the **EntropyEstimation** package.

## See also

[div_gen_simpson](https://ericmarcon.github.io/divent/dev/reference/div_gen_simpson.md)

\#' @references Grabchak M, Marcon E, Lang G, Zhang Z (2017). “The
Generalized Simpson's Entropy Is a Measure of Biodiversity.” *Plos One*,
**12**(3), e0173305.
[doi:10.1371/journal.pone.0173305](https://doi.org/10.1371/journal.pone.0173305)
.  
  
Zhang Z, Grabchak M (2016). “Entropic Representation and Estimation of
Diversity Indices.” *Journal of Nonparametric Statistics*, **28**(3),
563–575.
[doi:10.1080/10485252.2016.1190357](https://doi.org/10.1080/10485252.2016.1190357)
.  
  
Zhang Z, Zhou J (2010). “Re-Parameterization of Multinomial
Distributions and Diversity Indices.” *Journal of Statistical Planning
and Inference*, **140**(7), 1731–1738.
[doi:10.1016/j.jspi.2009.12.023](https://doi.org/10.1016/j.jspi.2009.12.023)
.

## Examples

``` r
# Entropy of each community
ent_gen_simpson(paracou_6_abd, k = 50)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 Zhang        50   0.472
#> 2 subplot_2   1.56 Zhang        50   0.528
#> 3 subplot_3   1.56 Zhang        50   0.503
#> 4 subplot_4   1.56 Zhang        50   0.498
# gamma entropy
ent_gen_simpson(paracou_6_abd, k = 50, gamma = TRUE)
#> # A tibble: 1 × 3
#>   estimator order entropy
#>   <chr>     <dbl>   <dbl>
#> 1 Zhang        50   0.515
```
