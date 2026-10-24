# Blinded Sample Size Re-Estimation for a Stratified Trial

Re-estimates the final group sizes of the strata of an ongoing two-arm
trial with a binary endpoint from the blinded data of an interim
analysis: the numbers of patients and the pooled numbers of responders
of the strata.

## Usage

``` r
StratBSSR(
  n1,
  n2,
  S,
  Delta.A,
  r,
  alpha,
  tar.power,
  Test,
  restricted = FALSE,
  N1 = NULL,
  N2 = NULL,
  frac = NULL,
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker"),
  effect = c("RD", "RR", "OR"),
  ss.method = c("exact", "standard"),
  ss.Test = Test,
  ss.alpha = alpha,
  ss.strata = c("stratified", "pooled"),
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50),
  N.min = NULL,
  N.max = NULL
)
```

## Arguments

- n1:

  Numbers of patients of group 1 observed in the strata at the interim
  analysis

- n2:

  Numbers of patients of group 2 observed in the strata at the interim
  analysis

- S:

  Numbers of responders observed in the strata at the interim analysis,
  pooled over both groups, one element per stratum

- Delta.A:

  Assumed treatment effect, on the scale given by `effect`, common to
  the strata and used to divide the blinded estimate of the pooled
  response probability of each stratum between the groups at the
  re-estimation

- r:

  Allocation ratio of group 1 to group 2 within the strata, a single
  value or one value per stratum. The blinded estimate of the pooled
  response probability of a stratum is divided between the groups with
  this ratio, also when the interim sizes are not in this ratio

- alpha:

  Level of significance of the final analysis

- tar.power:

  Target power of the re-estimation

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

- restricted:

  Logical. If `TRUE`, the final group sizes are at least the planned
  sizes. Default is `FALSE`

- N1:

  Planned sizes of group 1 in the strata, or `NULL` (default). Required
  when `restricted` is `TRUE`, when `frac` is `NULL`, and when the
  blinded data show no effect in the direction of the alternative

- N2:

  Planned sizes of group 2 in the strata, or `NULL` (default)

- frac:

  Fractions of the strata used to divide the re-estimated total sample
  size, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md).
  Default is `NULL`, which takes the fractions of the planned design

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

- effect:

  Scale of `Delta.A` and `Delta.T`. Options: `'RD'` (default) for the
  risk difference, `'RR'` for the risk ratio or `'OR'` for the odds
  ratio

- ss.method:

  Re-estimation rule. `'exact'` (default) re-estimates the sample size
  of
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)
  with `method = 'exact'`, and `'standard'` that of the normal
  approximation. With three or more strata only `'standard'` is
  available

- ss.Test:

  Test used in the re-estimation. Default is `Test`

- ss.alpha:

  Level used in the re-estimation. Default is `alpha`

- ss.strata:

  Blinded estimate of the re-estimation. `'stratified'` (default)
  estimates the pooled response probability of each stratum from its own
  interim patients, and `'pooled'` uses the proportion of responders
  among all the interim patients for every stratum, as a re-estimation
  that ignores the strata would. With `'pooled'` the final group sizes
  depend only on the total number of interim responders

- search:

  Search of the exact re-estimation, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)

- search.limit:

  Limit of the search `'stable'`, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)

- N.min:

  Smallest final total sample size, or `NULL` (default)

- N.max:

  Largest final total sample size, or `NULL` (default)

## Value

An object of class `bbssr_sbssr`, a data frame with one row per stratum
containing:

- Stratum:

  Number of the stratum

- n1, n2, S:

  Interim data

- hat.p:

  Blinded estimate of the pooled response probability

- hat.p1, hat.p2:

  Planning probabilities of the two groups, truncated to the unit
  interval

- N1.final, N2.final:

  Final group sizes

- n1.stage2, n2.stage2:

  Numbers of patients still to be enrolled

The attributes `n.raw` (unrounded total of the normal approximation),
`N.final` (final total sample size) and `Power` (exact power of the
final design at the planning probabilities) summarize the re-estimation,
and the attribute `search.limit` holds the limit of the nominal total
examined by `search = 'stable'` (otherwise `NA`).

## Details

The re-estimation is that of
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md):
in each stratum the blinded estimate of the pooled response probability
is divided between the groups with the assumed effect `Delta.A`, and the
total sample size of
[`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)
for these probabilities is divided between the strata. The final group
sizes are at least the interim sizes, at least the planned sizes when
`restricted = TRUE`, and their total lies between `N.min` and `N.max`.

## See also

[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
StratBSSR(n1 = c(15, 15), n2 = c(15, 15), S = c(12, 18), Delta.A = 0.2, r = 1,
          alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(30, 30),
          N2 = c(30, 30), ss.method = 'standard')
#> Blinded sample size re-estimation for a stratified trial
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Re-estimation    : normal approximation (unrestricted)
#>   Effect           : RD, assumed 0.2
#>   Final total      : 188, power 0.7985 at the planning probabilities
#> 
#>  Stratum n1 n2  S hat.p hat.p1 hat.p2 N1.final N2.final n1.stage2 n2.stage2
#>        1 15 15 12   0.4    0.5    0.3       47       47        32        32
#>        2 15 15 18   0.6    0.7    0.5       47       47        32        32
```
