# Plot phylo_divent Objects

Plot objects of class "phylo_divent" produced by
[as_phylo_divent](https://ericmarcon.github.io/divent/dev/reference/phylo_divent.md),
that are phylogenetic trees.

## Usage

``` r
# S3 method for class 'phylo_divent'
plot(x, ...)
```

## Arguments

- x:

  An object of class "phylo_divent".

- ...:

  Arguments passed to
  [stats::plot.dendrogram](https://rdrr.io/r/stats/dendrogram.html).

## Value

`NULL`. Called for side effects.

## Examples

``` r
# Paracou plot 6 species taxonomy
tree <- as_phylo_divent(paracou_6_taxo)
plot(tree, leaflab = "none")

```
