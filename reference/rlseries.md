# Log-Series Distribution

Random generation for the log-series distribution.

## Usage

``` r
rlseries(n, size, fisher_alpha, show_progress = TRUE, check_arguments = TRUE)
```

## Arguments

- n:

  the number of observations.

- size:

  The size of the distribution.

- fisher_alpha:

  Fisher's \\\alpha\\ in the log-series distribution.

- show_progress:

  if TRUE, a progress bar is shown during long computations.

- check_arguments:

  if `TRUE`, the function arguments are verified. Should be set to
  `FALSE` to save time when the arguments have been checked elsewhere.

## Value

A numeric vector with the random values drawn from the log-series
distribution.

## Details

Fast implementation of the random generation of a log-series
distribution (Fisher et al. 1943) .

The complete set of functions (including density, distribution function
and quantiles) can be found in package *sads* but this implementation of
the random generation is much faster.

If `size` is too large, i.e. `size` + 1 can't be distinguished from
`size` due to rounding, then an error is raised.

## References

Fisher RA, Corbet AS, Williams CB (1943). “The Relation between the
Number of Species and the Number of Individuals in a Random Sample of an
Animal Population.” *Journal of Animal Ecology*, **12**, 42–58.
[doi:10.2307/1411](https://doi.org/10.2307/1411) .

## Examples

``` r
# Generate a community made of 10000 individuals with alpha=40
size <- 1E4
fisher_alpha <- 40
species_number <- fisher_alpha * log(1 + size / fisher_alpha)
abundances <- rlseries(species_number, size = 1E5, fisher_alpha = 40)
# rcommunity() may be a better choice here
autoplot(rcommunity(1, size = 1E4, fisher_alpha = 40, distribution = "lseries"))
```
