# Similarity-Based Entropy of a Community

Estimate the entropy of species from abundance or probability data and a
similarity matrix between species. Several estimators are available to
deal with incomplete sampling. Bias correction requires the number of
individuals.

## Usage

``` r
ent_similarity(x, similarities, q = 1, ...)

# S3 method for class 'numeric'
ent_similarity(
  x,
  similarities = diag(length(x)),
  q = 1,
  estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang", "UnveilC", "UnveiliC",
    "naive"),
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  sample_coverage = NULL,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
ent_similarity(
  x,
  similarities = diag(sum(!colnames(x) %in% non_species_columns)),
  q = 1,
  estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang", "UnveilC", "UnveiliC",
    "naive"),
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
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
  If it is a numeric vector, then its length must equal the dimensions
  of the `similarities` matrix: species are assumed to be in the same
  order.

- similarities:

  a similarity matrix, that can be obtained by
  [fun_similarity](https://ericmarcon.github.io/divent/dev/reference/fun_similarity.md).
  Its default value is the identity matrix.

- q:

  a number: the order of diversity.

- ...:

  Unused.

- estimator:

  An estimator of entropy.

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

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/dev/reference/coverage.md).

- sample_coverage:

  the sample coverage of `x` calculated elsewhere. Used to calculate the
  gamma diversity of meta-communities, see details.

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

All species of the `species_distribution` must be found in the matrix of
`similarities` if it is named. If it is not or if `x` is numeric, its
size must equal the number of species. Then, the order of species is
assumed to be the same as that of the `species_distribution` or its
numeric equivalent.

Similarity-Based entropy can't be interpolated of extrapolated as of the
state of the art.

## Examples

``` r
# Similarity matrix
Z <- fun_similarity(paracou_6_fundist)
# Diversity of each community
ent_similarity(paracou_6_abd, similarities = Z, q = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       2   0.234
#> 2 subplot_2   1.56 UnveilJ       2   0.248
#> 3 subplot_3   1.56 UnveilJ       2   0.245
#> 4 subplot_4   1.56 UnveilJ       2   0.228
# gamma diversity
ent_similarity(paracou_6_abd, similarities = Z, q = 2, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity UnveilJ       2   0.237
```
