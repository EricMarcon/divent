# Allen et al.'s Phylogenetic Entropy of a Community

Estimate entropy (Allen et al. 2009) from abundance or probability data
and a phylogenetic or functional dendrogram.

## Usage

``` r
ent_allen(x, tree, ...)

# S3 method for class 'numeric'
ent_allen(
  x,
  tree,
  q = 1,
  normalize = TRUE,
  prune = FALSE,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
ent_allen(
  x,
  tree,
  q = 1,
  normalize = TRUE,
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

- q:

  a number: the order of diversity.

- normalize:

  if `TRUE`, phylogenetic is normalized: the height of the tree is set
  to 1.

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
entropy.

## Details

Estimators to deal with incomplete sampling are not implemented. Use
function
[ent_phylo](https://ericmarcon.github.io/divent/reference/ent_phylo.md)
with argument if they are needed.

The phylogenetic entropy is calculated following Allen et al. (2009) for
order \\q=1\\ and Leinster and Cobbold (2012) for other orders. The
result is identical to the total entropy calculated by
[ent_phylo](https://ericmarcon.github.io/divent/reference/ent_phylo.md).
It is much faster but no bias correction is available.

All species of the `species_distribution` must be found in the tips of
the `tree`.

## References

Allen B, Kon M, Bar-Yam Y (2009). “A New Phylogenetic Diversity Measure
Generalizing the Shannon Index and Its Application to Phyllostomid
Bats.” *American Naturalist*, **174**(2), 236–243.
[doi:10.1086/600101](https://doi.org/10.1086/600101) .  
  
Leinster T, Cobbold C (2012). “Measuring Diversity: The Importance of
Species Similarity.” *Ecology*, **93**(3), 477–489.
[doi:10.1890/10-2402.1](https://doi.org/10.1890/10-2402.1) .

## Examples

``` r
# entropy of each community
ent_allen(paracou_6_abd, tree = paracou_6_taxo)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 naive         1    3.60
#> 2 subplot_2   1.56 naive         1    3.83
#> 3 subplot_3   1.56 naive         1    3.74
#> 4 subplot_4   1.56 naive         1    3.63

# gamma entropy
ent_allen(paracou_6_abd, tree = paracou_6_taxo, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order entropy
#>   <chr>         <chr>     <dbl>   <dbl>
#> 1 Metacommunity naive         1    3.82
```
