# Fit a distribution

Fit a well-known distribution to a species distribution.

## Usage

``` r
fit_rac(x, ...)

# S3 method for class 'numeric'
fit_rac(
  x,
  distribution = c("lnorm", "lseries", "geom", "bstick"),
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
fit_rac(
  x,
  distribution = c("lnorm", "lseries", "geom", "bstick"),
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object

- ...:

  Unused.

- distribution:

  The distribution of species abundances. May be "lnorm" (log-normal),
  "lseries" (log-series), "geom" (geometric) or "bstick" (broken stick).

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A tibble with the sites and the estimated distribution parameters.

## Details

[abundances](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
can be used to fit rank-abundance curves (RAC) of classical
distributions:

- "lnorm" for log-normal (Preston 1948) .

- "lseries" for log-series (Fisher et al. 1943) .

- "geom" for geometric (Motomura 1932) .

- "bstick" for broken stick (MacArthur 1957) . It has no parameter, so
  the maximum abundance is returned.

## References

Fisher RA, Corbet AS, Williams CB (1943). “The Relation between the
Number of Species and the Number of Individuals in a Random Sample of an
Animal Population.” *Journal of Animal Ecology*, **12**, 42–58.
[doi:10.2307/1411](https://doi.org/10.2307/1411) .  
  
MacArthur RH (1957). “On the Relative Abundance of Bird Species.”
*Proceedings of the National Academy of Sciences of the United States of
America*, **43**(3), 293–295.
[doi:10.1073/pnas.43.3.293](https://doi.org/10.1073/pnas.43.3.293) ,
89566.  
  
Motomura I (1932). “On the statistical treatment of communities.”
*Zoological Magazine*, **44**, 379–383.  
  
Preston FW (1948). “The Commonness, and Rarity, of Species.” *Ecology*,
**29**(3), 254–283.
[doi:10.2307/1930989](https://doi.org/10.2307/1930989) .

## Examples

``` r
fit_rac(paracou_6_abd, distribution = "lnorm")
#> # A tibble: 4 × 3
#>   site         mu sigma
#>   <chr>     <dbl> <dbl>
#> 1 subplot_1 0.848 1.03 
#> 2 subplot_2 0.802 0.981
#> 3 subplot_3 0.848 1.00 
#> 4 subplot_4 0.823 0.979
```
