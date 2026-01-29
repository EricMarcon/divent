# Abundances of Communities

Utilities for community abundances (objects of class "abundances").

## Usage

``` r
abd_species(abundances, check_arguments = TRUE)

abd_sum(abundances, as_numeric = FALSE, check_arguments = TRUE)

prob_species(species_distribution, check_arguments = TRUE)
```

## Arguments

- abundances:

  an object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md).

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- species_distribution:

  an object of class
  [species_distribution](https://ericmarcon.github.io/divent/reference/species_distribution.md).

## Value

`abd_species()` returns a tibble containing the species abundance
columns only, to simplify numeric operations.

`prob_species()` returns the same tibble but values are probabilities.

`abd_sum()` returns the sample sizes of the communities in a numeric
vector.

## Examples

``` r
abd_species(paracou_6_abd)
#> # A tibble: 4 × 335
#>   Abarema_jupunba Abarema_mataybifolia Amaioua_guianensis Amanoa_congesta
#>             <int>                <int>              <int>           <int>
#> 1               2                    2                  1               1
#> 2               2                    0                  1               0
#> 3               2                    2                  0               0
#> 4               4                    0                  0               0
#> # ℹ 331 more variables: Amanoa_guianensis <int>, Ambelania_acida <int>,
#> #   Amphirrhox_longifolia <int>, Andira_coriacea <int>, Apeiba_glabra <int>,
#> #   Aspidosperma_album <int>, Aspidosperma_cruentum <int>,
#> #   Aspidosperma_excelsum <int>, Bocoa_prouacensis <int>,
#> #   Brosimum_guianense <int>, Brosimum_rubescens <int>, Brosimum_utile <int>,
#> #   Carapa_surinamensis <int>, Caryocar_glabrum <int>, Casearia_decandra <int>,
#> #   Casearia_javitensis <int>, Catostemma_fragrans <int>, …
prob_species(paracou_6_abd)
#> # A tibble: 4 × 335
#>   Abarema_jupunba Abarema_mataybifolia Amaioua_guianensis Amanoa_congesta
#>             <dbl>                <dbl>              <dbl>           <dbl>
#> 1         0.00212              0.00212            0.00106         0.00106
#> 2         0.00229              0                  0.00115         0      
#> 3         0.00215              0.00215            0               0      
#> 4         0.00501              0                  0               0      
#> # ℹ 331 more variables: Amanoa_guianensis <dbl>, Ambelania_acida <dbl>,
#> #   Amphirrhox_longifolia <dbl>, Andira_coriacea <dbl>, Apeiba_glabra <dbl>,
#> #   Aspidosperma_album <dbl>, Aspidosperma_cruentum <dbl>,
#> #   Aspidosperma_excelsum <dbl>, Bocoa_prouacensis <dbl>,
#> #   Brosimum_guianense <dbl>, Brosimum_rubescens <dbl>, Brosimum_utile <dbl>,
#> #   Carapa_surinamensis <dbl>, Caryocar_glabrum <dbl>, Casearia_decandra <dbl>,
#> #   Casearia_javitensis <dbl>, Catostemma_fragrans <dbl>, …
abd_sum(paracou_6_abd)
#> # A tibble: 4 × 3
#>   site      weight abundance
#>   <chr>      <dbl>     <dbl>
#> 1 subplot_1   1.56       942
#> 2 subplot_2   1.56       872
#> 3 subplot_3   1.56       929
#> 4 subplot_4   1.56       798
```
