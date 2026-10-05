# Plot Method for bbssr_grid Objects

Draws the power or the expected total sample size of the designs of
[`BinaryGridBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryGridBSSR.md)
against the true pooled response probability, with one colour per value
of a chosen column and, optionally, one panel per value of another.

## Usage

``` r
# S3 method for class 'bbssr_grid'
plot(
  x,
  what = c("power", "E.N"),
  colour.by = NULL,
  facet.by = NULL,
  main = NULL,
  sub = NULL,
  xlab = NULL,
  ylab = NULL,
  ref.line = NULL,
  legend.title = NULL,
  base_size = 11,
  ...
)
```

## Arguments

- x:

  An object of class `bbssr_grid`

- what:

  Quantity to draw. `'power'` (default) draws the power of the BSSR
  designs as solid lines and that of the fixed-sample designs as dashed
  lines. `'E.N'` draws the expected total sample size of the BSSR
  designs

- colour.by:

  Name of the column of `x` whose values are told apart by colour. The
  default of `NULL` uses the index `design`

- facet.by:

  Name of a column of `x` whose values are drawn in separate panels, or
  `NULL` (default) for a single panel

- main:

  Title. `NULL` keeps the default and `NA` removes it

- sub:

  Subtitle. `NULL` keeps the default and `NA` removes it

- xlab:

  Label of the horizontal axis. `NULL` keeps the default

- ylab:

  Label of the vertical axis. `NULL` keeps the default

- ref.line:

  Height of the dotted reference line for `what = 'power'`. `NULL` draws
  it at the target power when all designs share one, and `NA` removes it

- legend.title:

  Title of the colour legend. `NULL` uses `colour.by`

- base_size:

  Base font size passed to `theme_bw()`

- ...:

  Further arguments, ignored

## Value

A `ggplot` object

## Examples

``` r
design <- expand.grid(Test = c('Chisq', 'Fisher'), omega = c(0.3, 0.5),
                      stringsAsFactors = FALSE)
grid <- BinaryGridBSSR(design, p = c(0.3, 0.4, 0.5), Delta.A = 0.3, N1 = 30, N2 = 30,
                       r = 1, alpha = 0.025, tar.power = 0.8, ss.method = 'standard')
plot(grid, colour.by = 'Test', facet.by = 'omega')

plot(grid, what = 'E.N', colour.by = 'Test', facet.by = 'omega')

```
