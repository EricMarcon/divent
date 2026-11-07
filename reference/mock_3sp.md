# Mock data

A simple dataset to test diversity functions. It contains 3 species with
their abundances, their distance matrix and their phylogenetic tree.

## Usage

``` r
mock_3sp_abd

mock_3sp_dist

mock_3sp_tree
```

## Format

`mock_3sp_abd` is a vector, `mock_3sp_dist` a matrix and `mock_3sp_tree`
an object of class
[ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html).

An object of class `dist` of length 3.

An object of class `phylo` of length 4.

## Examples

``` r
mock_3sp_abd
#> A B C 
#> 3 2 1 
mock_3sp_dist
#>   A B
#> B 1  
#> C 2 2
if (require("ape")) {
  plot(mock_3sp_tree, direction = "downwards")
}
#> Loading required package: ape
#> 
#> Attaching package: ‘ape’
#> The following objects are masked from ‘package:spatstat.geom’:
#> 
#>     edges, rotate
axis(2)

```
