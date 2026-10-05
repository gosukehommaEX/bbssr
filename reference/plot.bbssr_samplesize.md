# Plot the Exact Power Curve around a Sample Size Solution

Recomputes the exact power over a range of sample sizes of group 2 and
marks the selected sample size and the target power.

## Usage

``` r
# S3 method for class 'bbssr_samplesize'
plot(
  x,
  N2.range = NULL,
  main = NULL,
  sub = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ref.line = NULL,
  ref.line.N2 = NULL,
  show.points = TRUE,
  colours = NULL,
  base_size = 11,
  ...
)
```

## Arguments

- x:

  An object of class `bbssr_samplesize` returned by
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)

- N2.range:

  Optional integer vector of sample sizes of group 2 at which the power
  is evaluated. By default the selected sample size plus or minus five
  is used

- main:

  Title of the plot. The default of `NULL` builds the title from the
  object and a single `NA` drops it

- sub:

  Subtitle of the plot, following the convention of `main`

- xlab:

  Label of the horizontal axis, following the convention of `main`

- ylab:

  Label of the vertical axis, following the convention of `main`

- ylim:

  Numeric vector of length two giving the range of the vertical axis, or
  `NULL` for a range chosen from the data

- ref.line:

  Numeric vector of heights at which a dashed horizontal line is drawn.
  The default of `NULL` places one line at the target power and a single
  `NA` draws no line

- ref.line.N2:

  Numeric vector of positions at which a dotted vertical line is drawn.
  The default of `NULL` places one line at the selected sample size and
  a single `NA` draws no line

- show.points:

  Logical. If `TRUE`, the default, a marker is drawn at every evaluated
  sample size. Set it to `FALSE` to draw the curve alone, which is
  easier to read over a wide range

- colours:

  Colour of the curve

- base_size:

  Base font size of the theme, in points

- ...:

  Further arguments, currently ignored

## Value

A `ggplot` object

## Details

The power is recomputed at every point of `N2.range`, so a wide range
combined with one of the unconditional tests can take a long time to
evaluate.

The vertical range is imposed with
[`coord_cartesian`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html),
so points outside `ylim` are hidden rather than removed from the data.

## Examples

``` r
# \donttest{
ss <- BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
                       tar.power = 0.8, Test = 'Chisq')
plot(ss)


# Enlarge the type, rename the axes and draw neither reference line
plot(ss, base_size = 14, ref.line = NA, ref.line.N2 = NA,
     xlab = 'Sample size of the control group', ylab = 'Exact power')

# }
```
