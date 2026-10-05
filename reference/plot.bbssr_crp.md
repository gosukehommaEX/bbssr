# Plot Method for bbssr_crp Objects

Draws the conditional rejection probabilities of a design with blinded
sample size re-estimation over the pooled responder counts of the two
stages, or the type I error rate given the pooled number of interim
responders.

## Usage

``` r
# S3 method for class 'bbssr_crp'
plot(
  x,
  what = c("CRP", "CRP.total", "by.s"),
  main = NULL,
  sub = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.title = NULL,
  colours = NULL,
  base_size = 11,
  ...
)
```

## Arguments

- x:

  An object of class `bbssr_crp`

- what:

  Quantity to draw. `'CRP'` (default) and `'CRP.total'` draw the
  corresponding column over the pairs `(s, s2)`, shaded from white to
  the first colour, and outline the outcomes at which it exceeds the
  nominal level in the second colour. `'by.s'` draws the type I error
  rate given `s` for each value of `theta`, which requires an object
  created with `theta`

- main:

  Title. `NULL` keeps the default and `NA` removes it

- sub:

  Subtitle. `NULL` keeps the default and `NA` removes it

- xlab:

  Label of the horizontal axis. `NULL` keeps the default

- ylab:

  Label of the vertical axis. `NULL` keeps the default

- legend.title:

  Title of the legend. `NULL` keeps the default

- colours:

  Two colours, for the shading of the probabilities and for the outline
  of the outcomes above the nominal level. They are not used for
  `what = 'by.s'`

- base_size:

  Base font size passed to `theme_bw()`

- ...:

  Further arguments, ignored

## Value

A `ggplot` object

## Examples

``` r
crp <- BinaryCondRejectBSSR(
  Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Fisher', ss.method = 'standard',
  theta = c(0.3, 0.5)
)
plot(crp)

plot(crp, what = 'by.s')

```
