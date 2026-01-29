# Paracou plot 6

A community assembly. It contains number of trees per species of the
plot \#6 of Paracou. The plot covers 6.25 ha of tropical rainforest,
divided into 4 equally-sized subplots.

## Usage

``` r
paracou_6_abd

paracou_6_wmppp
```

## Format

`paracou_6_abd` is an object of class
[abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md),
which is also a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html).
Each line of the tibble is a subplot. `paracou_6_wmppp` is a
[dbmss::wmppp](https://ericmarcon.github.io/dbmss/reference/wmppp.html)
object, i.e. a weighted, marked planar point pattern.

An object of class `wmppp` (inherits from `ppp`) of length 6.

## Source

Permanent data census of Paracou: <https://paracou.cirad.fr/>

## Details

In `paracou_6_abd` (a tibble), the "site" column contains the subplot
number, "weight" contains its area and all others columns contain a
species. Data are the number of trees above 10 cm diameter at breast
height (DBH).

In `paracou_6_wmppp` (a point pattern), the point type is tree species
and the point weight is their basal area, in square centimeters.

This dataset is from Paracou field station, French Guiana, managed by
[Cirad](https://www.cirad.fr).

## See also

[paracou_6_taxo](https://ericmarcon.github.io/divent/reference/paracou_6_taxo.md),
[paracou_6_fundist](https://ericmarcon.github.io/divent/reference/paracou_6_fundist.md)

## Examples

``` r
# Rank-abundance curve of the species of the whole plot
autoplot(metacommunity(paracou_6_abd))

```
