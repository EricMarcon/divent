# Phylogenetic Diversity Accumulation of a Community

Diversity and Entropy Accumulation Curves represent the accumulation of
entropy with respect to the sample size.

## Usage

``` r
accum_ent_phylo(x, ...)

# S3 method for class 'numeric'
accum_ent_phylo(
  x,
  tree,
  q = 0,
  normalize = TRUE,
  levels = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  n_simulations = 0,
  alpha = 0.05,
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'abundances'
accum_ent_phylo(
  x,
  tree,
  q = 0,
  normalize = TRUE,
  levels = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  gamma = FALSE,
  n_simulations = 0,
  alpha = 0.05,
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
)

accum_div_phylo(x, ...)

# S3 method for class 'numeric'
accum_div_phylo(
  x,
  tree,
  q = 0,
  normalize = TRUE,
  levels = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  n_simulations = 0,
  alpha = 0.05,
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
)

# S3 method for class 'abundances'
accum_div_phylo(
  x,
  tree,
  q = 0,
  normalize = TRUE,
  levels = NULL,
  probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
  unveiling = c("geometric", "uniform", "none"),
  richness_estimator = c("rarefy", "jackknife", "iChao1", "Chao1", "naive"),
  jack_alpha = 0.05,
  jack_max = 10,
  coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
  gamma = FALSE,
  n_simulations = 0,
  alpha = 0.05,
  show_progress = TRUE,
  ...,
  check_arguments = TRUE
)
```

## Arguments

- x:

  An object, that may be a numeric vector containing abundances or
  probabilities, or an object of class
  [abundances](https://ericmarcon.github.io/divent/reference/species_distribution.md)
  or
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md).

- ...:

  Unused.

- tree:

  an ultrametric, phylogenetic tree. May be an object of class
  [phylo_divent](https://ericmarcon.github.io/divent/reference/phylo_divent.md),
  [ape::phylo](https://rdrr.io/pkg/ape/man/read.tree.html),
  [ade4::phylog](https://adeverse.github.io/ade4/reference/phylog.html)
  or [stats::hclust](https://rdrr.io/r/stats/hclust.html).

- q:

  a number: the order of diversity.

- normalize:

  if `TRUE`, phylogenetic is normalized: the height of the tree is set
  to 1.

- levels:

  The levels, i.e. the sample sizes of interpolation or extrapolation: a
  vector of integer values.

- probability_estimator:

  a string containing one of the possible estimators of the probability
  distribution (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only for extrapolation.

- unveiling:

  a string containing one of the possible unveiling methods to estimate
  the probabilities of the unobserved species (see
  [probabilities](https://ericmarcon.github.io/divent/reference/probabilities.md)).
  Used only for extrapolation.

- richness_estimator:

  an estimator of richness to evaluate the total number of species, see
  [div_richness](https://ericmarcon.github.io/divent/reference/div_richness.md).
  used for interpolation and extrapolation.

- jack_alpha:

  the risk level, 5% by default, used to optimize the jackknife order.

- jack_max:

  the highest jackknife order allowed. Default is 10.

- coverage_estimator:

  an estimator of sample coverage used by
  [coverage](https://ericmarcon.github.io/divent/reference/coverage.md).

- n_simulations:

  the number of simulations used to estimate the confidence envelope.

- alpha:

  the risk level, 5% by default.

- show_progress:

  if TRUE, a progress bar is shown during long computations.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

- gamma:

  if `TRUE`, \\\gamma\\ diversity, i.e. diversity of the metacommunity,
  is computed.

## Value

A tibble with the site names, the estimators used and the accumulated
entropy or diversity at each level of sampling effort.

## Details

`accum_ent_phylo()` or `accum_div_phylo()` estimate the phylogenetic
diversity or entropy accumulation curve of a distribution. See
[ent_tsallis](https://ericmarcon.github.io/divent/reference/ent_tsallis.md)
for details about the computation of entropy at each level of
interpolation and extrapolation.

In accumulation curves, extrapolation if done by estimating the
asymptotic distribution of the community and estimating entropy at
different levels by interpolation.

Interpolation and extrapolation of integer orders of diversity are from
Chao et al. (2014) . The asymptotic richness is adjusted so that the
extrapolated part of the accumulation joins the observed value at the
sample size.

"accumulation" objects can be plotted. They generalize the classical
Species Accumulation Curves (SAC) which are diversity accumulation of
order \\q=0\\.

## References

Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM
(2014). “Rarefaction and Extrapolation with Hill Numbers: A Framework
for Sampling and Estimation in Species Diversity Studies.” *Ecological
Monographs*, **84**(1), 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) .

## Examples

``` r
# Richness accumulation up to the sample size.
# 100 simulations only to save time.
autoplot(
  accum_div_phylo(mock_3sp_abd, tree = mock_3sp_tree, n_simulations = 100)
)

```
