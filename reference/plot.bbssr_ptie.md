# Plot Method for bbssr_ptie Objects

Draws the type I error rate of a matched-pairs design with blinded
sample size re-estimation and of the fixed design against the
probability of a discordant pair on the boundary of the null hypothesis,
with a reference line at the nominal level. For a test of equivalence
the two boundaries are drawn in separate panels.

## Usage

``` r
# S3 method for class 'bbssr_ptie'
plot(
  x,
  main = NULL,
  sub = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ref.line = NULL,
  legend.title = NULL,
  legend.labels = NULL,
  colours = NULL,
  base_size = 11,
  ...
)
```

## Arguments

- x:

  An object of class `bbssr_ptie`

- main:

  Title. `NULL` keeps the default and `NA` removes it

- sub:

  Subtitle. `NULL` keeps the default and `NA` removes it

- xlab:

  Label of the horizontal axis. `NULL` keeps the default

- ylab:

  Label of the vertical axis. `NULL` keeps the default

- ylim:

  Range of the vertical axis, applied through `coord_cartesian()`

- ref.line:

  Height of the dashed reference line. `NULL` draws it at the nominal
  level and `NA` removes it

- legend.title:

  Title of the legend. `NULL` keeps the default

- legend.labels:

  Labels of the two designs. `NULL` keeps the defaults

- colours:

  Two colours for the BSSR and the fixed design

- base_size:

  Base font size passed to `theme_bw()`

- ...:

  Further arguments, ignored

## Value

A `ggplot` object

## Examples

``` r
tie <- PairedTypeIErrorBSSR(Delta.A = 0.2, N = 40, n.interim = 20, alpha = 0.05,
                            tar.power = 0.8, Test = 'McNemar',
                            alternative = 'two.sided', ss.method = 'standard',
                            psi = seq(0.1, 0.9, by = 0.05))
plot(tie)

```
