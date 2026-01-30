# Rao's Quadratic Entropy of a Community

Estimate the quadratic entropy (Rao 1982) of species from abundance or
probability data. An estimator (Lande 1996) is available to deal with
incomplete sampling.

## Usage

``` r
ent_rao(x, ...)

# S3 method for class 'numeric'
ent_rao(
  x,
  distances = NULL,
  tree = NULL,
  normalize = TRUE,
  estimator = c("Lande", "naive"),
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
ent_rao(
  x,
  distances = NULL,
  tree = NULL,
  normalize = TRUE,
  estimator = c("Lande", "naive"),
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
  [abundances](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md).

- ...:

  Unused.

- distances:

  a distance matrix or an object of class
  [stats::dist](https://rdrr.io/r/stats/dist.html).

- tree:

  an ultrametric, phylogenetic tree. May be an object of class
  [phylo_divent](https://ericmarcon.github.io/divent/dev/reference/phylo_divent.md),
  [ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html),
  [ade4::phylog](https://adeverse.github.io/ade4/reference/phylog.html)
  or [stats::hclust](https://rdrr.io/r/stats/hclust.html).

- normalize:

  if `TRUE`, phylogenetic is normalized: the height of the tree is set
  to 1.

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

Rao's entropy is phylogenetic or similarity-based entropy of order 2.
[ent_phylo](https://ericmarcon.github.io/divent/dev/reference/ent_phylo.md)
and
[ent_similarity](https://ericmarcon.github.io/divent/dev/reference/ent_similarity.md)
with argument `q = 2` provide more estimators and allow estimating
entropy at a chosen level.

All species of the `species_distribution` must be found in the matrix of
`distances` if it is named. If it is not or if `x` is numeric, its size
must equal the number of species. Then, the order of species is assumed
to be the same as that of the `species_distribution` or its numeric
equivalent.

## References

Lande R (1996). “Statistics and Partitioning of Species Diversity, and
Similarity among Multiple Communities.” *Oikos*, **76**(1), 5–13.
[doi:10.2307/3545743](https://doi.org/10.2307/3545743) .  
  
Rao CR (1982). “Diversity and Dissimilarity Coefficients: A Unified
Approach.” *Theoretical Population Biology*, **21**, 24–43.
[doi:10.1016/0040-5809(82)90004-1](https://doi.org/10.1016/0040-5809%2882%2990004-1)
.

## Examples

``` r
# Entropy of each community
ent_rao(paracou_6_abd, tree = paracou_6_taxo)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 Lande         2   0.970
#> 2 subplot_2   1.56 Lande         2   0.977
#> 3 subplot_3   1.56 Lande         2   0.973
#> 4 subplot_4   1.56 Lande         2   0.973
# Similar to (but estimators are not the same)
ent_phylo(paracou_6_abd, tree = paracou_6_taxo, q = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator     q entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       2   0.943
#> 2 subplot_2   1.56 UnveilJ       2   0.953
#> 3 subplot_3   1.56 UnveilJ       2   0.951
#> 4 subplot_4   1.56 UnveilJ       2   0.939

# Functional entropy
ent_rao(paracou_6_abd, distances = paracou_6_fundist)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 Lande         2   0.365
#> 2 subplot_2   1.56 Lande         2   0.393
#> 3 subplot_3   1.56 Lande         2   0.383
#> 4 subplot_4   1.56 Lande         2   0.365

# gamma entropy
ent_rao(paracou_6_abd, tree = paracou_6_taxo, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity Lande         2   0.976
```
