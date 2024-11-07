# Diversity measures based on Entropy<img src="man/figures/logo.png" align="right" alt="" width="120" />

![stability-wip](https://img.shields.io/badge/stability-work_in_progress-lightgrey.svg)
![R-CMD-check](https://github.com/EricMarcon/divent/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/github/EricMarcon/divent/branch/master/graphs/badge.svg)](https://app.codecov.io/github/EricMarcon/divent)
[![CodeFactor](https://www.codefactor.io/repository/github/ericmarcon/divent/badge/master)](https://www.codefactor.io/repository/github/ericmarcon/divent/overview/master)
[![CRAN version](https://www.r-pkg.org/badges/version/divent)](https://CRAN.r-project.org/package=divent)
[![](https://cranlogs.r-pkg.org/badges/grand-total/divent)](https://CRAN.R-project.org/package=divent)
[![](https://cranlogs.r-pkg.org/badges/divent)](https://CRAN.R-project.org/package=divent)

**divent** is an R package that provides functions to estimate alpha, beta and gamma diversity of communities, including phylogenetic and functional diversity.

It is a reboot of the package **entropart** to make it tidy, easier to use and optimize the code that has been added along years of research.

## Details

In the divent package, individuals of different *species* are counted in several *communities* which may (or not) be agregated to define a *metacommunity*. 
In the metacommunity, the probability to find a species in the weighted average of probabilities in communities. 
This is a naming convention, which may correspond to plots in a forest inventory or any data organized the same way.

Basic functions allow computing diversity of a community.
Calculate entropies such as Shannon's, Simpson's, or Hurlbert's and explicit diversity (i.e. effective number of species), aka Hill numbers.
By default, the best available estimator of diversity will be used, according to the data.
  
Communities can be simulated and plotted.
  
Phylogenetic entropy and diversity can be calculated if a phylogenetic (or functional), ultrametric tree is provided, with the state-of-the-art estimation-bias correction. 
Similarity-based diversity is calculated, based on a similarity matrix.

The diversity of a metacommunity, i.e. $\gamma$ diversity, can be partitioned into $\alpha$ (that of communities) and $\beta$ diversities.

# Vignettes

A quick [introduction](https://ericmarcon.github.io/divent/articles/divent.html) is in `vignette("divent")`.

A full documentation is available online, in the "Articles" section of the web site of the vignette.

The [documentation of the development version](https://EricMarcon.github.io/divent/dev/) is also available.


## Reference

Marcon, E. and Herault, B. (2015). entropart: An R Package to Measure and Partition Diversity.
*Journal of Statistical Software*. 67(8): 1-26.
