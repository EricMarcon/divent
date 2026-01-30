# Hill number of a Community

Estimate the diversity sensu stricto, i.e. the Hill (1973) number of
species from abundance or probability data.

## Usage

``` r
div_hill(x, q = 1, ...)

# S3 method for class 'numeric'
div_hill(
  x,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q_threshold = 10,
  sample_coverage = NULL,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'species_distribution'
div_hill(
  x,
  q = 1,
  estimator = c("UnveilJ", "ChaoJost", "ChaoShen", "GenCov", "Grassberger", "Marcon",
    "UnveilC", "UnveiliC", "ZhangGrabchak", "naive", "Bonachela", "Holste"),
  level = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  q_threshold = 10,
  gamma = FALSE,
  as_numeric = FALSE,
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object, that may be a numeric vector containing abundances or
  probabilities, or an object of class
  [abundances](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md).

- q:

  a number: the order of diversity.

- ...:

  Unused.

- estimator:

  an estimator of asymptotic diversity.

- level:

  the level of interpolation or extrapolation. It may be a sample size
  (an integer) or a sample coverage (a number between 0 and 1). If not
  `NULL`, the asymptotic `estimator` is ignored.

- probability_estimator:

  a string containing one of the possible estimators of the probability
  distribution (see
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md)).
  Used only for extrapolation.

- unveiling:

  a string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species (see
  [probabilities](https://ericmarcon.github.io/divent/dev/reference/probabilities.md)).
  Used only for extrapolation.

- richness_estimator:

  an estimator of richness to evaluate the total number of species, see
  [div_richness](https://ericmarcon.github.io/divent/dev/reference/div_richness.md).
  used for interpolation and extrapolation.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/dev/reference/coverage.md).

- q_threshold:

  the value of `q` above which diversity is computed directly with the
  naive estimator \\(\sum{p_s^q}^{\frac{1}{(1-q)}}\\, without computing
  entropy. When `q` is great, the exponential of entropy goes to
  \\0^{\frac{1}{(1-q)}}\\, causing rounding errors while the naive
  estimator of diversity is less and less biased.

- sample_coverage:

  the sample coverage of `x` calculated elsewhere. Used to calculate the
  gamma diversity of meta-communities, see details.

- as_numeric:

  if `TRUE`, a number or a numeric vector is returned rather than a
  tibble.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- gamma:

  if `TRUE`, \\\gamma\\ diversity, i.e. diversity of the metacommunity,
  is computed.

## Value

A tibble with the site names, the estimators used and the estimated
diversity.

## Details

Several estimators are available to deal with incomplete sampling.

Bias correction requires the number of individuals.

Estimation techniques are from Chao and Shen (2003) , Grassberger (1988)
,Holste et al. (1998) , Bonachela et al. (2008) , Marcon et al. (2014)
which is actually the max value of "ChaoShen" and "Grassberger", Zhang
and Grabchak (2014) , Chao et al. (2015) , Chao and Jost (2015) and
Marcon (2015) .

The `ChaoJost` estimator (Chao et al. 2013; Chao and Jost 2015) contains
an unbiased part concerning observed species, equal to that of Zhang and
Grabchak (2014) , and a (biased) estimator of the remaining bias based
on the estimation of the species-accumulation curve. It is very
efficient but slow if the number of individuals is more than a few
hundreds.

The unveiled estimators rely on Chao et al. (2015) , completed by Marcon
(2015) . The actual probabilities of observed species are estimated and
completed by a geometric distribution of the probabilities of unobserved
species. The number of unobserved species is estimated by the Chao1
estimator (`UnveilC`), following Chao et al. (2015) , or by the iChao1
(`UnveiliC`) or the jackknife (`UnveilJ`). The `UnveilJ` estimator often
has a lower bias but a greater variance (Marcon 2015) . It is a good
first choice thanks to the versatility of the jackknife estimator of
richness.

Estimators by Bonachela et al. (2008) and Holste et al. (1998) are
rarely used.

To estimate \\\gamma\\ diversity, the size of a metacommunity (see
[metacommunity](https://ericmarcon.github.io/divent/dev/reference/metacommunity.md))
is unknown so it has to be set according to a rule which does not ensure
that its abundances are integer values. Then, classical bias-correction
methods do not apply. Providing the `sample_coverage` argument allows
applying the `ChaoShen` and `Grassberger` estimators to estimate quite
well the entropy.

Diversity can be estimated at a specified level of interpolation or
extrapolation, either a chosen sample size or sample coverage (Chao et
al. 2014) , rather than its asymptotic value. See
[accum_hill](https://ericmarcon.github.io/divent/dev/reference/accum_hill.md)
for details.

## References

Bonachela JA, Hinrichsen H, Muñoz MA (2008). “Entropy Estimates of Small
Data Sets.” *Journal of Physics A: Mathematical and Theoretical*,
**41**(202001), 1–9.
[doi:10.1088/1751-8113/41/20/202001](https://doi.org/10.1088/1751-8113/41/20/202001)
.  
  
Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .  
  
Chao A, Hsieh TC, Chazdon RL, Colwell RK, Gotelli NJ (2015). “Unveiling
the Species-Rank Abundance Distribution by Generalizing Good-Turing
Sample Coverage Theory.” *Ecology*, **96**(5), 1189–1201.
[doi:10.1890/14-0550.1](https://doi.org/10.1890/14-0550.1) .  
  
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
[doi:10.1111/2041-210x.12108](https://doi.org/10.1111/2041-210x.12108)
.  
  
Grassberger P (1988). “Finite Sample Corrections to Entropy and
Dimension Estimates.” *Physics Letters A*, **128**(6-7), 369–373.
[doi:10.1016/0375-9601(88)90193-4](https://doi.org/10.1016/0375-9601%2888%2990193-4)
.  
  
Hill MO (1973). “Diversity and Evenness: A Unifying Notation and Its
Consequences.” *Ecology*, **54**(2), 427–432.
[doi:10.2307/1934352](https://doi.org/10.2307/1934352) .  
  
Holste D, Große I, Herzel H (1998). “Bayes' Estimators of Generalized
Entropies.” *Journal of Physics A: Mathematical and General*,
**31**(11), 2551–2566.  
  
Marcon E (2015). “Practical Estimation of Diversity from Abundance
Data.” *HAL*, **01212435**(version 2).  
  
Marcon E, Scotti I, Hérault B, Rossi V, Lang G (2014). “Generalization
of the Partitioning of Shannon Diversity.” *Plos One*, **9**(3), e90289.
[doi:10.1371/journal.pone.0090289](https://doi.org/10.1371/journal.pone.0090289)
.  
  
Zhang Z, Grabchak M (2014). “Nonparametric Estimation of
Kullback-Leibler Divergence.” *Neural computation*, **26**(11),
2570–2593.
[doi:10.1162/NECO_a_00646](https://doi.org/10.1162/NECO_a_00646) ,
25058703.

## Examples

``` r
# Diversity of each community
div_hill(paracou_6_abd, q = 2)
#> # A tibble: 4 × 5
#>   site      weight estimator order diversity
#>   <chr>      <dbl> <chr>     <dbl>     <dbl>
#> 1 subplot_1   1.56 UnveilJ       2      42.3
#> 2 subplot_2   1.56 UnveilJ       2      44.6
#> 3 subplot_3   1.56 UnveilJ       2      48.8
#> 4 subplot_4   1.56 UnveilJ       2      36.0
# gamma diversity
div_hill(paracou_6_abd, q = 2, gamma = TRUE)
#> # A tibble: 1 × 4
#>   site          estimator order diversity
#>   <chr>         <chr>     <dbl>     <dbl>
#> 1 Metacommunity UnveilJ       2      46.5

# At 80% coverage
div_hill(paracou_6_abd, q = 2, level = 0.8)
#> # A tibble: 4 × 6
#>   site      weight estimator order level diversity
#>   <chr>      <dbl> <chr>     <dbl> <dbl>     <dbl>
#> 1 subplot_1   1.56 Chao2014      2   304      37.2
#> 2 subplot_2   1.56 Chao2014      2   347      39.6
#> 3 subplot_3   1.56 Chao2014      2   333      42.7
#> 4 subplot_4   1.56 Chao2014      2   303      32.3
```
