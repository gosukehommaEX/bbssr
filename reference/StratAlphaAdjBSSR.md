# Adjusted Significance Level of Blinded Sample Size Re-Estimation for a Stratified Design

Lowers the nominal significance level of a two-stratum design with
blinded sample size re-estimation until the largest type I error rate
over the common response probabilities of the strata does not exceed the
target level `alpha`, and returns the adjusted level together with that
of the fixed design with the planned group sizes.

## Usage

``` r
StratAlphaAdjBSSR(
  Delta.A,
  N1,
  N2,
  omega = NULL,
  r,
  alpha,
  tar.power,
  Test,
  restricted = FALSE,
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
  N.max = NULL,
  n.interim = NULL,
  frac = NULL,
  theta = seq(0, 1, by = 0.05),
  maximize = c("certified", "grid"),
  adjust = c("test", "both"),
  tol = 1e-08,
  step = 1e-05,
  alloc = c("fixed", "stratum", "group", "both"),
  block = 4,
  ss.frac = c("observed", "planned"),
  nseq = 100
)
```

## Arguments

- Delta.A:

  Assumed treatment effect, on the scale given by `effect`, common to
  the strata and used to divide the blinded estimate of the pooled
  response probability of each stratum between the groups at the
  re-estimation

- N1:

  Planned sizes of group 1 in the strata, which are those of the fixed
  design

- N2:

  Planned sizes of group 2 in the strata

- omega:

  Fraction of the planned sizes at the interim analysis: in stratum `j`
  the interim size of group 2 is `ceiling(omega N2[j])` and that of
  group 1 is `ceiling(r[j] n2[j])`. Supply either `omega` or `n.interim`

- r:

  Allocation ratio of group 1 to group 2 within the strata, a single
  value or one value per stratum. The blinded estimate of the pooled
  response probability of a stratum is divided between the groups with
  this ratio, also when the interim sizes are not in this ratio

- alpha:

  Target level of significance, which the type I error rate must not
  exceed

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

  Level of significance used for the re-estimation when
  `adjust = 'test'`. Default is `alpha`. It is ignored when
  `adjust = 'both'`

- ss.strata:

  Blinded estimate of the re-estimation. `'stratified'` (default)
  estimates the pooled response probability of each stratum from its own
  interim patients, and `'pooled'` uses the proportion of responders
  among all the interim patients for every stratum, as a re-estimation
  that ignores the strata would. With `'pooled'` the final group sizes
  depend only on the total number of interim responders

- search:

  Search of the exact re-estimation, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md);
  it does not go beyond the largest size that `N.max` allows

- search.limit:

  Limit of the search `'stable'`, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)

- N.min:

  Smallest final total sample size, or `NULL` (default)

- N.max:

  Largest final total sample size, or `NULL` (default)

- n.interim:

  Interim group sizes, as a matrix with two rows (group 1 and group 2)
  and one column per stratum. Supply either `omega` or `n.interim`

- frac:

  Fractions of the strata used to divide the re-estimated total sample
  size between the strata, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md).
  Default is `NULL`, which takes the fractions of the planned design,
  `(N1 + N2) / sum(N1 + N2)`

- theta:

  Common response probabilities of the strata under the null hypothesis.
  A vector gives, with two strata, the grid of all pairs of its values,
  and with three or more strata the points at which every stratum has
  the same value. A matrix with one column per stratum gives the points
  directly. Default is `seq(0, 1, by = 0.05)`

- maximize:

  How the largest type I error rate at a level is located, as in
  [`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md).
  With `'certified'` (default) the adjusted level controls the type I
  error rate over the whole rectangle spanned by `theta`. With `'grid'`
  the control is assessed at the points of `theta`

- adjust:

  Which part of the design uses the adjusted level. `'test'` (default)
  applies it to the final analysis only, so the re-estimation keeps
  using `ss.alpha`. `'both'` applies it to the re-estimation as well

- tol:

  Tolerance of the bisection, relative to `alpha`, a single value in (0,
  1). Default is `1e-8`

- step:

  Decrement of the level in the search used when `adjust = 'both'`.
  Default is `1e-5`

- alloc:

  Allocation of the patients to the strata and the groups. `'fixed'`
  (default) uses the planned and interim group sizes. The others make
  the sizes random, see Details: `'stratum'` (stratum sizes random,
  permuted blocks within the strata) and `'both'` (stratum and group
  sizes random), with any number of strata, and `'group'` (stratum sizes
  in the fixed fractions, group sizes random), with two strata

- block:

  Size of the permuted blocks of `alloc = 'stratum'`. Each block holds
  `block r / (1 + r)` patients of group 1, which must be a whole number.
  Default is 4

- ss.frac:

  Fractions of the strata used by the re-estimation with a random
  allocation: `'observed'` (default), those of the interim patients, or
  `'planned'`, those given by `frac`

- nseq:

  Number of allocation sequences drawn with a random allocation and two
  strata. Default is 100

## Value

An object of class `bbssr_alphaadj`, a data frame with one row for the
BSSR design and one for the fixed design containing:

- Design:

  `'BSSR'` or `'TRAD'`

- alpha:

  Target level

- max.TIE:

  Largest type I error rate at the nominal level `alpha`

- alpha.adj:

  Adjusted nominal level. It is `alpha` when the nominal level already
  controls the type I error rate, and 0, with a warning, when no level
  examined by the bisection does. The search of `adjust = 'both'` for
  the BSSR design stops with an error instead when it reaches 0

- max.TIE.adj:

  Largest type I error rate at the adjusted level

- theta.S1.adj, theta.S2.adj:

  Common response probabilities of the strata at which `max.TIE.adj`
  occurs, `NA` when `alpha.adj` is 0

- max.TIE.bound:

  Upper bound of the type I error rate at the nominal level with
  `maximize = 'certified'`, otherwise `NA`

- max.TIE.adj.bound:

  Upper bound of the type I error rate at the adjusted level with
  `maximize = 'certified'`, otherwise `NA`

## Details

The searches are those of
[`BinaryAlphaAdjBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md).
With `adjust = 'test'` the final group sizes do not depend on the level,
the largest type I error rate is a non-decreasing step function of the
level, and a bisection locates the level above which it exceeds `alpha`,
to within `tol * alpha`. The reported level is the largest value with at
most six significant digits that gives the same critical values. With
`adjust = 'both'` the level is lowered from `alpha` in steps of `step`
until the type I error rate is controlled. Each level requires the
coefficients of the type I error rate of
[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md),
so the search takes about as long as that many evaluations of the type I
error rate.

With a random allocation (`alloc` other than `'fixed'`) the type I error
rate of the design with re-estimation at each level is that of
[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md),
the average over the same `nseq` allocation sequences at every level, so
each level takes about `nseq` times as long as with the fixed
allocation. The adjusted level controls this average, whose Monte Carlo
error it does not account for. The rate of the fixed design with the
planned total and the same allocation is exact. Both levels are reported
as the search found them, without the rounding to six significant
digits, since the critical values differ between the allocation
sequences and between the stratum and group sizes.

## See also

[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
# \donttest{
StratAlphaAdjBSSR(Delta.A = 0.25, N1 = c(15, 15), N2 = c(15, 15), omega = 0.5, r = 1,
                  alpha = 0.025, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
                  theta = seq(0.05, 0.95, by = 0.05))
#> Adjusted significance level for blinded sample size re-estimation
#> 
#>   Test            : CMH
#>   Alternative     : greater
#>   Adjusted part   : final analysis only
#>   Maximum         : certified over (theta.S1, theta.S2) in [0.05, 0.95] x [0.05, 0.95]
#>   Target level    : 0.025
#> 
#>        Design   max.TIE alpha.adj max.TIE.adj
#>          BSSR 0.0268795 0.0242744   0.0247366
#>  Fixed sample 0.0277886 0.0209052   0.0235488
# }
```
