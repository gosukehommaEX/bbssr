# Plot the Power of a Stratified BSSR Design against the Fixed Design

Draws the power of the design with blinded sample size re-estimation of
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md)
and that of the fixed design against a chosen column of the scenarios,
as
[`plot.bbssr_powerbssr`](https://gosukehommaex.github.io/bbssr/reference/plot.bbssr_powerbssr.md)
does for two groups.

## Usage

``` r
# S3 method for class 'bbssr_spowerbssr'
plot(x, x.by = NULL, xlab = NULL, ...)
```

## Arguments

- x:

  An object of class `bbssr_spowerbssr` returned by
  [`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md)

- x.by:

  Name of the numeric column of `x` on the horizontal axis, or
  `'scenario'` for the rows of `p` that the scenarios come from. The
  default of `NULL` takes the first of the true pooled response
  probabilities `p.S1, p.S2, ...` that varies between the scenarios, or
  the scenario when none varies

- xlab:

  Label of the horizontal axis. `NULL` makes it from `x.by` and a single
  `NA` drops it

- ...:

  Further arguments of
  [`plot.bbssr_powerbssr`](https://gosukehommaex.github.io/bbssr/reference/plot.bbssr_powerbssr.md):
  `main`, `sub`, `ylab`, `ylim`, `ref.line`, `legend.title`,
  `legend.labels`, `show.points`, `colours` and `base_size`

## Value

A `ggplot` object

## Details

The Monte Carlo standard errors of a simulated or averaged power (three
or more strata, or a random allocation) are not drawn. With
`Delta.T = 0` in
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md)
the curves are rejection probabilities under the null hypothesis, and
`ref.line` can be set to the level of significance.

## Examples

``` r
res <- StratPowerBSSR(p = cbind(0.4, seq(0.3, 0.7, by = 0.1)), Delta.A = 0.3,
                      Delta.T = 0.3, N1 = c(10, 10), N2 = c(10, 10), omega = 0.5,
                      r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                      ss.method = 'standard')
plot(res)

```
