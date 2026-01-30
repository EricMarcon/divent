# Class phylo_divent

Methods for dendrograms of class "phylo_divent".

## Usage

``` r
as_phylo_divent(tree)

is_phylo_divent(x)
```

## Arguments

- tree:

  an ultrametric, phylogenetic tree. May be an object of class
  phylo_divent,
  [ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html),
  [ade4::phylog](https://adeverse.github.io/ade4/reference/phylog.html)
  or [stats::hclust](https://rdrr.io/r/stats/hclust.html).

- x:

  An object of class "phylo_divent".

## Value

`as_phylo_divent` returns a phylogenetic tree that is an object of class
"phylo_divent".

## Details

`as_phylo_divent` calculates cuts and intervals of a phylogenetic tree
and makes it available both in
[stats::hclust](https://rdrr.io/r/stats/hclust.html) and
[ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html) formats. The
conversion preprocesses the tree: it calculates cuts so that the tree
can be reused efficiently by phylodiversity functions.

## Examples

``` r
# Paracou plot 6 species taxonomy
tree <- as_phylo_divent(mock_3sp_tree)
plot(tree)

```
