# Probabilities of Species

Estimate actual probabilities of species from a sample

## Usage

``` r
probabilities(x, ...)

# S3 method for class 'numeric'
probabilities(
  x,
  estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
  unveiling = c("none", "uniform", "geometric"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q = 0,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'abundances'
probabilities(
  x,
  estimator = c("naive", "Chao2013", "Chao2015", "ChaoShen"),
  unveiling = c("none", "uniform", "geometric"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "rarefy", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q = 0,
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object. It may be:

  - a numeric vector containing abundances. It may be named to track
    species names.

  - an object of class
    [species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md).

- ...:

  Unused.

- estimator:

  One of the estimators of a probability distribution: "naive" (the
  default value), or "Chao2013", "Chao2015", "ChaoShen" to estimate the
  probabilities of the observed species in the asymptotic distribution.

- unveiling:

  One of the possible unveiling methods to estimate the probabilities of
  the unobserved species: "none" (default, no species is added),
  "uniform" (all unobserved species have the same probability) or
  "geometric" (the unobserved species distribution is geometric).

- richness_estimator:

  An estimator of richness to evaluate the total number of species.
  "jackknife" is the default value. An alternative is "rarefy" to
  estimate the number of species such that the entropy of the asymptotic
  distribution rarefied to the observed sample size equals the actual
  entropy of the data.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/dev/reference/coverage.md).

- q:

  The order of diversity. Default is 0 for richness. Used only to
  estimate asymptotic probability distributions when argument
  `richness_estimator` is "rarefy". Then, the number of unobserved
  species is fitted so that the entropy of order q of the asymptotic
  probability distribution at the observed sample size equals the actual
  entropy of the data.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

An object of class "probabilities", which is a
[species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
or a numeric vector with argument `as_numeric = TRUE`.

## Details

`probabilities()` estimates a probability distribution from a sample. If
the `estimator` is not "naive", the observed abundance distribution is
used to estimate the actual probability distribution. The list of
species is changed: zero-abundance species are cleared, and some
unobserved species are added. First, observed species probabilities are
estimated following Chao and Shen (2003) , i.e. input probabilities are
multiplied by the sample coverage, or according to more sophisticated
models: Chao et al. (2013) , a single-parameter model, or Chao and Jost
(2015) , a two-parameter model. The total probability of observed
species equals the sample coverage. Then, the distribution of unobserved
species can be unveiled: their number is estimated according to the
`richness_estimator`. The coverage deficit (1 minus the sample coverage)
is shared by the unobserved species equally (`unveiling` = "uniform",
(Chao et al. 2013) ) or according to a geometric distribution
(`unveiling` = "geometric", (Chao and Jost 2015) ).

## References

Chao A, Jost L (2015). “Estimating Diversity and Entropy Profiles via
Discovery Rates of New Species.” *Methods in Ecology and Evolution*,
**6**(8), 873–882.
[doi:10.1111/2041-210X.12349](https://doi.org/10.1111/2041-210X.12349)
.  
  
Chao A, Shen T (2003). “Nonparametric Estimation of Shannon's Index of
Diversity When There Are Unseen Species in Sample.” *Environmental and
Ecological Statistics*, **10**(4), 429–443.
[doi:10.1023/A:1026096204727](https://doi.org/10.1023/A%3A1026096204727)
.  
  
Chao A, Wang Y, Jost L (2013). “Entropy and the Species Accumulation
Curve: A Novel Entropy Estimator via Discovery Rates of New Species.”
*Methods in Ecology and Evolution*, **4**(11), 1091–1100.
[doi:10.1111/2041-210x.12108](https://doi.org/10.1111/2041-210x.12108) .

## Examples

``` r
# Just transform abundances into probabilities
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
# Estimate the distribution of probabilities from observed abundances (unveiled probabilities)
prob_unv <- probabilities(
  paracou_6_abd,
  estimator = "Chao2015",
  unveiling = "geometric",
  richness_estimator = "jackknife"
)
# Estimate entropy from the unveiled probabilities
ent_shannon(prob_unv)
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1  1.00  naive         1    4.57
#> 2 subplot_2  1.00  naive         1    4.73
#> 3 subplot_3  1.000 naive         1    4.65
#> 4 subplot_4  1     naive         1    4.55
# Identical to
ent_shannon(paracou_6_abd, estimator = "UnveilJ")
#> # A tibble: 4 × 5
#>   site      weight estimator order entropy
#>   <chr>      <dbl> <chr>     <dbl>   <dbl>
#> 1 subplot_1   1.56 UnveilJ       1    4.57
#> 2 subplot_2   1.56 UnveilJ       1    4.73
#> 3 subplot_3   1.56 UnveilJ       1    4.65
#> 4 subplot_4   1.56 UnveilJ       1    4.55
```
