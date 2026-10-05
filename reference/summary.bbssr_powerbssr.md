# Summary Method for bbssr_powerbssr Objects

Summarizes the distribution of the final total sample size of a design
with blinded sample size re-estimation for every true pooled response
probability.

## Usage

``` r
# S3 method for class 'bbssr_powerbssr'
summary(object, probs = c(0.25, 0.5, 0.75), ...)
```

## Arguments

- object:

  An object of class `bbssr_powerbssr`

- probs:

  Probabilities of the reported quantiles. Default is
  `c(0.25, 0.5, 0.75)`

- ...:

  Further arguments, ignored

## Value

A data frame with one row per element of `p` containing `p1`, `p2`, `p`,
the expected final total sample size `E.N`, its standard deviation
`SD.N`, the quantiles of the final total sample size in columns named
after `probs` (for example `N.q50` for the median), and the probability
`P.N.max` that the final total reaches its largest attainable value

## Details

The quantile for the probability `q` is the smallest final total sample
size whose cumulative probability is at least `q`.

## Examples

``` r
res <- BinaryPowerBSSR(
  p = c(0.3, 0.45),
  Delta.A = 0.3, Delta.T = 0.3,
  N1 = 10, N2 = 10, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
)
summary(res)
#>     p1   p2    p      E.N     SD.N N.q25 N.q50 N.q75    P.N.max
#> 1 0.45 0.15 0.30 65.30594 15.41083    48    68    80 0.09908785
#> 2 0.60 0.30 0.45 75.57870 12.17796    68    80    80 0.24579081
```
