# Non-Species columns

The column names that are not those of species in a
[species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md).

## Usage

``` r
non_species_columns
```

## Format

A character vector.

## Details

Objects of classes
[abundances](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
and
[probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md),
that are also of class
[species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md),
have columns named after the species they contain. Some columns are
reserved to describe the plots and their diversity.

## Examples

``` r
# Default values
non_species_columns
#> [1] "site"      "weight"    "cut"       "interval"  "q"         "entropy"  
#> [7] "diversity" "abundance" "comments" 
```
