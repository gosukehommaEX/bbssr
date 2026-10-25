# Summary Method for bbssr_spowerbssr Objects

Summarizes the distribution of the final total sample size of a
stratified design with blinded sample size re-estimation in each
scenario of
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md).

## Usage

``` r
# S3 method for class 'bbssr_spowerbssr'
summary(object, probs = c(0.25, 0.5, 0.75), ...)
```

## Arguments

- object:

  An object of class `bbssr_spowerbssr`

- probs:

  Probabilities of the reported quantiles. Default is
  `c(0.25, 0.5, 0.75)`

- ...:

  Further arguments, ignored

## Value

A data frame with one row per scenario containing the true pooled
response probabilities `p.S1, p.S2, ...` of the strata, those of the two
groups `p1.S1, ...` and `p2.S1, ...`, the expected value `E.N` and the
standard deviation `SD.N` of the final total sample size, its quantiles
(`N.q25`, `N.q50` and `N.q75` by default, the smallest total whose
distribution function reaches the probability), and the probability
`P.N.max` that it is the largest total the design can reach. With two
strata and the fixed allocation this is the largest final total of the
re-estimation over all interim outcomes. Otherwise it is `N.max`, and
`P.N.max` is `NA` when `N.max` was not given.

## Details

The distribution is the attribute `N.dist` of `object`: exact with two
strata and the fixed allocation, the average over the allocation
sequences with two strata and a random allocation, and that of the
simulated trials with three or more strata.

## Examples

``` r
res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10),
                      N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05,
                      tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
summary(res)
#>   p.S1 p.S2 p1.S1 p1.S2 p2.S1 p2.S2      E.N     SD.N N.q25 N.q50 N.q75
#> 1  0.4  0.6  0.55  0.75  0.25  0.45 59.74021 7.128156    56    60    64
#>     P.N.max
#> 1 0.1985219
```
