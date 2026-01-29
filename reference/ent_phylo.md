# Phylogenetic Entropy of a Community

Estimate the entropy of species from abundance or probability data and a
phylogenetic tree. Several estimators are available to deal with
incomplete sampling.

## Usage

``` r
ent_phylo(x, tree, q = 1, ...)

# S3 method for class 'numeric'
ent_phylo(
  x,
  tree,
  q = 1,
  normalize = TRUE,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
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
ent_phylo(
  x,
  tree,
  q = 1,
  normalize = TRUE,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
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

  An object, that may be a named numeric vector (names are species
  names) containing abundances or probabilities, or an object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md).

- tree:

  an ultrametric, phylogenetic tree. May be an object of class
  [phylo_divent](https://ericmarcon.github.io/divent/reference/phylo_divent.md),
  [ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html),
  [ade4::phylog](https://adeverse.github.io/ade4/reference/phylog.html)
  or [stats::hclust](https://rdrr.io/r/stats/hclust.html).

- q:

  a number: the order of diversity.

- ...:

  Unused.

- normalize:

  if `TRUE`, phylogenetic is normalized: the height of the tree is set
  to 1.

- estimator:

  An estimator of entropy.

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

Bias correction requires the number of individuals. See
[div_hill](https://ericmarcon.github.io/divent/reference/div_hill.md)
for estimators.

Entropy can be estimated at a specified level of interpolation or
extrapolation, either a chosen sample size or sample coverage (Chao et
al. 2014) , rather than its asymptotic value. See
[accum_tsallis](https://ericmarcon.github.io/divent/reference/accum_hill.md)
for details.

## References

Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .

## Examples

``` r
# Entropy of each community
ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator     q entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       2   0.943
#> 2 subplot_2   1.56 UnveilJ       2   0.953
#> 3 subplot_3   1.56 UnveilJ       2   0.951
#> 4 subplot_4   1.56 UnveilJ       2   0.939
# Gamma entropy
ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2, gamma = TRUE)
#> ! The estimator can't be applied to non-integer values.
#> ! The estimator can't be applied to non-integer values.
#> # A tibble: 1 × 4
#>   site          estimator     q entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity UnveilJ       2   0.949

# At 80% coverage
ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2, level = 0.8)
#> # A tibble: 4 × 5
#>   site      weight estimator     q entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       2   0.931
#> 2 subplot_2   1.56 UnveilJ       2   0.944
#> 3 subplot_3   1.56 UnveilJ       2   0.940
#> 4 subplot_4   1.56 UnveilJ       2   0.929

```
