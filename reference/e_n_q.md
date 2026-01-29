# Grassberger's expectation of n^q

Expected value of \\n^q\\ when \\n\\ follows a Poisson distribution of
parameter \\n\\.

## Usage

``` r
e_n_q(n, q)
```

## Arguments

- n:

  A positive integer vector.

- q:

  A positive number.

## Value

A vector of the same length as n containing the transformed values.

## Details

The expectation of \\n^q\\ when \\n\\ follows a Poisson distribution was
derived by Grassberger (1988) .

It is computed using the [beta](https://rdrr.io/r/base/Special.html)
function. Its value is 0 for \\n-q+1\<0\\, and close to 0 when \\n=q\\,
which is not a correct estimate: it should not be used when \\q \> n\\.

## References

Grassberger P (1988). “Finite Sample Corrections to Entropy and
Dimension Estimates.” *Physics Letters A*, **128**(6-7), 369–373.
[doi:10.1016/0375-9601(88)90193-4](https://doi.org/10.1016/0375-9601%2888%2990193-4)
.

## Examples

``` r
n <- 10
q <- 2
# Compare
e_n_q(n, q)
#> [1] 90
# with (empirical estimation)
mean(rpois(1000, lambda = n)^q)
#> [1] 107.415
# and (naive estimation)
n^q
#> [1] 100
```
