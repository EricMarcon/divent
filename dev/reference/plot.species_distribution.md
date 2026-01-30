# Plot Profile Objects

Plot objects of class "species_distribution" produced by
[species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md)
and similar functions.

## Usage

``` r
# S3 method for class 'species_distribution'
plot(
  x,
  type = c("RAC", "Metacommunity"),
  ...,
  fit_rac = FALSE,
  distribution = c("lnorm", "lseries", "geom", "bstick"),
  ylog = "y",
  main = NULL,
  xlab = "Rank",
  ylab = NULL,
  palette = "Set1"
)

# S3 method for class 'species_distribution'
autoplot(
  object,
  ...,
  fit_rac = FALSE,
  distribution = c("lnorm", "lseries", "geom", "bstick"),
  ylog = TRUE,
  main = NULL,
  xlab = "Rank",
  ylab = NULL,
  pch = 19,
  cex = 1.5
)
```

## Arguments

- x:

  An object.

- type:

  The type of plot. "RAC" (Rank-abundance curve, or Whittaker plot) or
  "Metacommunity" to represent species abundances of each community
  along with those of the metacommunity.

- ...:

  Additional arguments to be passed to
  [plot](https://rdrr.io/r/graphics/plot.default.html). Unused
  elsewhere.

- fit_rac:

  If `TRUE`, estimate a theoretical distribution and fit the data with
  it. RAC plot only.

- distribution:

  The distribution of species abundances. May be "lnorm" (log-normal),
  "lseries" (log-series), "geom" (geometric) or "bstick" (broken stick).
  RAC plot only.

- ylog:

  If `TRUE`, the Y-axis is in log-scale. RAC plot only.

- main:

  The title of the plot.

- xlab:

  The label of the X-axis. RAC plot only.

- ylab:

  The label of the Y-axis.

- palette:

  The name of a color palette, recognized by
  [RColorBrewer::brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).
  RAC plot only.

- object:

  An object of class
  [species_distribution](https://ericmarcon.github.io/divent/dev/reference/species_distribution.md).

- pch:

  The plotting characters. See
  [graphics::points](https://rdrr.io/r/graphics/points.html).

- cex:

  The character expansion (size) of the points. See
  [graphics::points](https://rdrr.io/r/graphics/points.html).

## Value

`NULL`. Called for side effects.
