# Exact Power of a Stratified Test

Calculates the exact power of a test of
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md)
for a two-arm trial with a binary endpoint and fixed group sizes in each
stratum.

## Usage

``` r
StratPower(
  p1,
  p2,
  N1,
  N2,
  alpha,
  Test,
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker")
)
```

## Arguments

- p1:

  Response probabilities of group 1 in the strata: a vector with one
  element per stratum, or a matrix with one row per scenario and one
  column per stratum

- p2:

  Response probabilities of group 2 in the strata, of the same form as
  `p1`

- N1:

  Sizes of group 1 in the strata, one element per stratum

- N2:

  Sizes of group 2 in the strata, one element per stratum

- alpha:

  Level of significance, for the alternative specified

- Test:

  Statistical test. Options: `'CMH'` (Cochran-Mantel-Haenszel test with
  the hypergeometric variance, without continuity correction),
  `'CMH-cc'` (the same with the continuity correction of
  [`mantelhaen.test`](https://rdrr.io/r/stats/mantelhaen.test.html)),
  `'Cochran'` (Cochran's test with the binomial variance, the test of
  Woolson et al. (1986)), `'Cochran-cc'` (the same with the continuity
  correction, which rejects the same outcomes as that of Nam (1992) at
  one-sided levels up to 0.078), `'Strat-Fisher'` (the stratified exact
  conditional test of Jung (2014)), or one of the unstratified tests of
  the table of totals: `'Unstrat-Chisq'` (Pearson's chi-squared test,
  the comparison of two proportions with the pooled variance),
  `'Unstrat-Chisq-cc'` (the same with the continuity correction) and
  `'Unstrat-Fisher'` (Fisher's exact test). See Details

- alternative:

  Direction of the alternative hypothesis. Options: `'greater'`
  (default), under which group 1 responds more often than group 2,
  `'less'` or `'two.sided'`

- tsmethod:

  Convention of the two-sided p-value of `'Strat-Fisher'` and
  `'Unstrat-Fisher'`, as for the Fisher exact test of
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md).
  Options: `'minlike'` (default), `'central'` or `'blaker'`. Ignored by
  the other tests and for one-sided alternatives

## Value

An object of class `bbssr_spower`, a data frame with one row per
scenario containing the response probabilities `p1.S1, p1.S2, ...` and
`p2.S1, p2.S2, ...` and the exact power `Power`.

## Details

The numbers of responders of the two groups in each stratum are
independent binomial variables. The power is the probability of the
rejection region of
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md),
computed exactly by enumerating the numbers of responders of the strata
and carrying the distribution of the partial sum of `S` over the strata.
The response probabilities may differ between the strata in any way, so
the power is available also when the odds ratio or the risk difference
is not common to the strata.

## See also

[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md),
[`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), N1 = c(40, 40), N2 = c(40, 40),
           alpha = 0.025, Test = 'CMH')
#> Exact power of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Group sizes      : N1 = (40, 40), N2 = (40, 40)
#>   Alpha            : 0.025
#> 
#>  p1.S1 p1.S2 p2.S1 p2.S2  Power
#>    0.5   0.7   0.3   0.5 0.7322

# Several scenarios, and the stratified exact test
StratPower(p1 = rbind(c(0.5, 0.7), c(0.45, 0.65)), p2 = rbind(c(0.3, 0.5), c(0.3, 0.5)),
           N1 = c(40, 40), N2 = c(40, 40), alpha = 0.025, Test = 'Strat-Fisher')
#> Exact power of a stratified test
#> 
#>   Test             : Strat-Fisher
#>   Alternative      : greater
#>   Group sizes      : N1 = (40, 40), N2 = (40, 40)
#>   Alpha            : 0.025
#> 
#>  p1.S1 p1.S2 p2.S1 p2.S2  Power
#>   0.50  0.70   0.3   0.5 0.6792
#>   0.45  0.65   0.3   0.5 0.4249
```
