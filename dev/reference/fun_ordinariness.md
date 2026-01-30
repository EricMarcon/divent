# Functional ordinariness of a community

The ordinariness of a species is the average similarity of its
individuals with others (Leinster and Cobbold 2012) .

## Usage

``` r
fun_ordinariness(
  species_distribution,
  similarities = diag(sum(!colnames(species_distribution) %in% non_species_columns)),
  as_numeric = FALSE,
  check_arguments = TRUE
)
```

## Arguments

- species_distribution:

  an object of class
  [species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md).

- similarities:

  a similarity matrix, that can be obtained by
  [fun_similarity](https://ericmarcon.github.io/divent/dev/reference/fun_similarity.md).
  Its default value is the identity matrix.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A tibble with the ordinariness of each species, or a matrix if argument
`as_numeric` is `TRUE`.

## Details

All species of the `species_distribution` must be found in the matrix of
`similarities` if it is named. If it is not, its size must equal the
number of species. Then, the order of species is assumed to be the same
as that of the `species_distribution`.

## References

Leinster T, Cobbold C (2012). “Measuring Diversity: The Importance of
Species Similarity.” *Ecology*, **93**(3), 477–489.
[doi:10.1890/10-2402.1](https://doi.org/10.1890/10-2402.1) .

## Examples

``` r
fun_ordinariness(paracou_6_abd, fun_similarity(paracou_6_fundist))
#> # A tibble: 4 × 337
#>   site      weight Abarema_jupunba Abarema_mataybifolia Amaioua_guianensis
#>   <chr>      <dbl>           <dbl>                <dbl>              <dbl>
#> 1 subplot_1   1.56           0.794                0.801              0.702
#> 2 subplot_2   1.56           0.790                0.785              0.694
#> 3 subplot_3   1.56           0.790                0.793              0.701
#> 4 subplot_4   1.56           0.797                0.797              0.695
#> # ℹ 332 more variables: Amanoa_congesta <dbl>, Amanoa_guianensis <dbl>,
#> #   Ambelania_acida <dbl>, Amphirrhox_longifolia <dbl>, Andira_coriacea <dbl>,
#> #   Apeiba_glabra <dbl>, Aspidosperma_album <dbl>, Aspidosperma_cruentum <dbl>,
#> #   Aspidosperma_excelsum <dbl>, Bocoa_prouacensis <dbl>,
#> #   Brosimum_guianense <dbl>, Brosimum_rubescens <dbl>, Brosimum_utile <dbl>,
#> #   Carapa_surinamensis <dbl>, Caryocar_glabrum <dbl>, Casearia_decandra <dbl>,
#> #   Casearia_javitensis <dbl>, Catostemma_fragrans <dbl>, …

# Compare with probabilities
probabilities(paracou_6_abd)
#> # A tibble: 4 × 337
#>   site      weight Abarema_jupunba Abarema_mataybifolia Amaioua_guianensis
#>   <chr>      <dbl>           <dbl>                <dbl>              <dbl>
#> 1 subplot_1  1             0.00212              0.00212            0.00106
#> 2 subplot_2  1.000         0.00229              0                  0.00115
#> 3 subplot_3  1.00          0.00215              0.00215            0      
#> 4 subplot_4  1.000         0.00501              0                  0      
#> # ℹ 332 more variables: Amanoa_congesta <dbl>, Amanoa_guianensis <dbl>,
#> #   Ambelania_acida <dbl>, Amphirrhox_longifolia <dbl>, Andira_coriacea <dbl>,
#> #   Apeiba_glabra <dbl>, Aspidosperma_album <dbl>, Aspidosperma_cruentum <dbl>,
#> #   Aspidosperma_excelsum <dbl>, Bocoa_prouacensis <dbl>,
#> #   Brosimum_guianense <dbl>, Brosimum_rubescens <dbl>, Brosimum_utile <dbl>,
#> #   Carapa_surinamensis <dbl>, Caryocar_glabrum <dbl>, Casearia_decandra <dbl>,
#> #   Casearia_javitensis <dbl>, Catostemma_fragrans <dbl>, …
# Decrease similarities so that ordinariness is close to probability
fun_ordinariness(paracou_6_abd, fun_similarity(paracou_6_fundist, rate = 100))
#> # A tibble: 4 × 337
#>   site      weight Abarema_jupunba Abarema_mataybifolia Amaioua_guianensis
#>   <chr>      <dbl>           <dbl>                <dbl>              <dbl>
#> 1 subplot_1   1.56         0.00217            0.00221          0.00110    
#> 2 subplot_2   1.56         0.00233            0.0000526        0.00116    
#> 3 subplot_3   1.56         0.00220            0.00225          0.00000939 
#> 4 subplot_4   1.56         0.00508            0.0000470        0.000000101
#> # ℹ 332 more variables: Amanoa_congesta <dbl>, Amanoa_guianensis <dbl>,
#> #   Ambelania_acida <dbl>, Amphirrhox_longifolia <dbl>, Andira_coriacea <dbl>,
#> #   Apeiba_glabra <dbl>, Aspidosperma_album <dbl>, Aspidosperma_cruentum <dbl>,
#> #   Aspidosperma_excelsum <dbl>, Bocoa_prouacensis <dbl>,
#> #   Brosimum_guianense <dbl>, Brosimum_rubescens <dbl>, Brosimum_utile <dbl>,
#> #   Carapa_surinamensis <dbl>, Caryocar_glabrum <dbl>, Casearia_decandra <dbl>,
#> #   Casearia_javitensis <dbl>, Catostemma_fragrans <dbl>, …
```
