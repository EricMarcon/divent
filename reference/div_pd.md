# Faith's Phylogenetic Diversity of a Community

Estimate PD (Faith 1992) or FD (Petchey and Gaston 2002) from abundance
or probability data and a phylogenetic or functional dendrogram.

## Usage

``` r
div_pd(x, tree, ...)

# S3 method for class 'numeric'
div_pd(x, tree, prune = FALSE, as_numeric = FALSE, ..., check_arguments = TRUE)

# S3 method for class 'species_distribution'
div_pd(
  x,
  tree,
  prune = FALSE,
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

- ...:

  Unused.

- prune:

  What to do when some species are in the tree but not in `x`? If
  `TRUE`, the tree is pruned to keep species of `x` only. The height of
  the tree may be changed if a pruned branch is related to the root. If
  `FALSE` (default), the length of branches of missing species is not
  summed but the height of the tree is never changed.

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
diversity.

## Details

Estimators to deal with incomplete sampling are not implemented. Use
function
[div_hill](https://ericmarcon.github.io/divent/reference/div_hill.md)
with argument `q = 0` if they are needed.

PD and FD are defined as the total length of the branches of the tree.

All species of the `species_distribution` must be found in the tips of
the `tree`.

## References

Faith DP (1992). “Conservation Evaluation and Phylogenetic Diversity.”
*Biological Conservation*, **61**(1), 1–10.
[doi:10.1016/0006-3207(92)91201-3](https://doi.org/10.1016/0006-3207%2892%2991201-3)
.  
  
Petchey OL, Gaston KJ (2002). “Functional Diversity (FD), Species
Richness and Community Composition.” *Ecology Letters*, **5**, 402–411.
[doi:10.1046/j.1461-0248.2002.00339.x](https://doi.org/10.1046/j.1461-0248.2002.00339.x)
.

## Examples

``` r
# diversity of each community
div_pd(paracou_6_abd, tree = paracou_6_taxo)
#> # A tibble: 4 × 5
#>   site      weight estimator order diversity
#>   <chr>      <dbl> <chr>     <dbl>     <dbl>
#> 1 subplot_1   1.56 naive         0       337
#> 2 subplot_2   1.56 naive         0       373
#> 3 subplot_3   1.56 naive         0       350
#> 4 subplot_4   1.56 naive         0       339

# gamma diversity
div_pd(paracou_6_abd, tree = paracou_6_taxo, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order diversity
#>   <chr>         <chr>     <dbl>     <dbl>
#> 1 Metacommunity naive         0       555
```
