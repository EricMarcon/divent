# alpha-shape calculation

Calculate a window containing all points of a point pattern. The window
is not convex but as close as possible to the points.

## Usage

``` r
alphahull(X, alpha = NULL)
```

## Arguments

- X:

  a planar point pattern
  ([spatstat.geom::ppp.object](https://rdrr.io/pkg/spatstat.geom/man/ppp.object.html)).

- alpha:

  a smoothing parameter to delimit concave polygons.

## Value

A window, i.e. a
[spatstat.geom::owin.object](https://rdrr.io/pkg/spatstat.geom/man/owin.object.html).

## Details

The typical use of this function is to define a narrow window around a
point pattern that has been created with a default, rectangle window.

The window is built by the
[`alphahull::ashape()`](https://rdrr.io/pkg/alphahull/man/ashape.html)
function and then transformed into a
[spatstat.geom::owin.object](https://rdrr.io/pkg/spatstat.geom/man/owin.object.html).
The `alpha` parameter determines the smallest size of zones excluded
from the window. If it is not specified, a first attempt is 1/256 of the
diameter of the existing window of `X`. If the shape cannot be
calculated, `alpha` is doubled and a new attempt is made.

## See also

[spatstat.geom::convexhull](https://rdrr.io/pkg/spatstat.geom/man/convexhull.html)

## Examples

``` r
# Simulate a point pattern
if (require(spatstat.random)) {
  X <- rpoispp(50)
  plot(X)
  # Calculate its border
  X$window <- alphahull(X)
  plot(X)
}
#> Loading required package: spatstat.random
#> Loading required package: spatstat.data
#> Loading required package: spatstat.univar
#> spatstat.univar 3.1-6
#> Loading required package: spatstat.geom
#> spatstat.geom 3.7-0
#> spatstat.random 3.4-4

#> Error in as_igraph_vs(graph, to): Invalid vertex names
```
