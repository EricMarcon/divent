# Aggregate communities into a metacommunity

Abundances of communities are summed according to their weights to
obtain the abundances of the metacommunity.

## Usage

``` r
metacommunity(x, name = "metacommunity", ...)

# S3 method for class 'matrix'
metacommunity(
  x,
  name = "metacommunity",
  weights = rep(1, nrow(x)),
  as_numeric = TRUE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'abundances'
metacommunity(
  x,
  name = "metacommunity",
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  that contains several communities or a matrix of abundances with
  communities in rows and species in columns.

- name:

  The name of the metacommunity

- ...:

  Unused.

- weights:

  the weights of the sites of the species distributions.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

An object of class
[abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md)
with a single row or a named vector if `as_numeric = TRUE`.

## Details

The total abundance of the metacommunity is by design equal to the sum
of community abundances so that the information used by diversity
estimators. A consequence is that equal weights lead to a metacommunity
whose species abundances are the sum of community species abundances.

If community weights are not equal then the metacommunity abundances are
in general not integer. Most diversity estimators can't be applied to
non-integer abundances but the knowledge of the sample coverage of each
community allow "ChaoShen" and "Grassberger" estimators.

## Examples

``` r
metacommunity(paracou_6_abd)
#> # A tibble: 1 × 337
#>   site          weight Abarema_jupunba Abarema_mataybifolia Amaioua_guianensis
#>   <chr>          <dbl>           <dbl>                <dbl>              <dbl>
#> 1 metacommunity   6.25              10                    4                  2
#> # ℹ 332 more variables: Amanoa_congesta <dbl>, Amanoa_guianensis <dbl>,
#> #   Ambelania_acida <dbl>, Amphirrhox_longifolia <dbl>, Andira_coriacea <dbl>,
#> #   Apeiba_glabra <dbl>, Aspidosperma_album <dbl>, Aspidosperma_cruentum <dbl>,
#> #   Aspidosperma_excelsum <dbl>, Bocoa_prouacensis <dbl>,
#> #   Brosimum_guianense <dbl>, Brosimum_rubescens <dbl>, Brosimum_utile <dbl>,
#> #   Carapa_surinamensis <dbl>, Caryocar_glabrum <dbl>, Casearia_decandra <dbl>,
#> #   Casearia_javitensis <dbl>, Catostemma_fragrans <dbl>, …
```
