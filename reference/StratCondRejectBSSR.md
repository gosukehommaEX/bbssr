# Conditional Rejection Probabilities of Blinded Sample Size Re-Estimation for a Stratified Design

Computes, under the null hypothesis of equal response probabilities in
every stratum, the probability that a trial with two strata and blinded
sample size re-estimation rejects the null hypothesis given the pooled
numbers of responders of the strata at the interim analysis and in the
second stage. These conditional rejection probabilities do not depend on
the response probabilities of the strata, and the type I error rate is
their average over the binomial distributions of the pooled counts. With
`theta`, the type I error rate is also returned for each interim
outcome. This is the counterpart of
[`BinaryCondRejectBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryCondRejectBSSR.md)
for a stratified design.

## Usage

``` r
StratCondRejectBSSR(
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
  theta = NULL,
  full = FALSE
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

  Optional common response probabilities of the strata at which the type
  I error rate is decomposed by the interim outcome: a vector, whose
  values form a grid of pairs, or a matrix with two columns that gives
  the points directly. Default is `NULL`, which skips the decomposition.
  The time of the decomposition grows with the number of points

- full:

  Logical. If `TRUE`, the attribute `full` lists the conditional
  rejection probability of every interim and second-stage outcome.
  Default is `FALSE`. A listing of more than 1e6 rows is refused

## Value

An object of class `bbssr_scrp`, a data frame with one row for every
pair of interim pooled numbers of responders of the strata, the first
stratum varying fastest, containing:

- s1, s2:

  Pooled numbers of responders of the two strata at the interim analysis

- N1.S1, N1.S2, N2.S1, N2.S2:

  Final group sizes reached from `(s1, s2)`

- CRP.max, CRP.min:

  Largest and smallest rejection probability over the pooled numbers of
  responders `(w1, w2)` of the strata in the second stage

- exceed:

  Number of pairs `(w1, w2)` at which the rejection probability exceeds
  `alpha`

- CRP.total.max:

  Largest conditional size of the final test given the final numbers of
  responders of the strata, see Details

The attribute `max` is a data frame with the largest values of `CRP` and
`CRP.total` over all outcomes, and the first outcome, in the order of
the rows and then of `(w1, w2)` with `w1` varying fastest, at which each
is attained up to rounding error. The attribute `exceed` gives the
numbers of outcomes `(s1, s2, w1, w2)` at which each of the two exceeds
`alpha`, and the attribute `reestimation` holds the final group sizes
for every interim outcome, as in
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md).
When `theta` is supplied, the attribute `by.s` is a data frame with one
row for each point and interim outcome, containing the probability
`prob.s` of `(s1, s2)`, the type I error rate `TIE.s` given `(s1, s2)`
and their product `contribution`, and the attribute `TIE` gives the type
I error rate at each point. With `full = TRUE` the attribute `full` is a
data frame with one row for every outcome `(s1, s2, w1, w2)`, ordered by
`s2`, `s1`, `w2` and `w1`, with the columns `CRP` and `CRP.total`.

## Details

Let `n1[j]` and `n2[j]` be the interim group sizes of stratum `j`, and
`m1[j]` and `m2[j]` the second-stage group sizes reached from the
interim outcome `s = (s1, s2)`. Under the null hypothesis both groups of
stratum `j` respond with the probability `theta_j`. Given `s_j`, the
number of interim responders of group 1 in stratum `j` is
hypergeometric, being the number of group 1 patients among `s_j`
responders drawn from the `n1[j] + n2[j]` interim patients of the
stratum. Given the pooled number `w_j` of second-stage responders of the
stratum, the number of second-stage responders of group 1 is an
independent hypergeometric count of the same kind. `CRP` sums the
rejection region of the final design over these four distributions and
does not depend on `theta`. The type I error rate is the sum over `s`
and `w = (w1, w2)` of
`b(s1; n1[1] + n2[1], theta_1) b(s2; n1[2] + n2[2], theta_2) b(w1; m1[1] + m2[1], theta_1) b(w2; m1[2] + m2[2], theta_2) CRP(s, w)`,
where `b` is the binomial probability.

`CRP.total` is the conditional size of the final test given the final
numbers of responders `z_j = s_j + w_j` of the strata, with the number
of responders of group 1 of each stratum hypergeometric among all its
final patients. The sum of the contributions at a point equals the type
I error rate of
[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md)
there.

Only designs with two strata are covered, as for the exact computations
of
[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md).

## See also

[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md),
[`BinaryCondRejectBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryCondRejectBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
crp <- StratCondRejectBSSR(Delta.A = 0.3, N1 = c(8, 8), N2 = c(8, 8), omega = 0.5,
                           r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                           ss.method = 'standard', theta = c(0.3, 0.5))
print(crp)
#> Conditional rejection probabilities of blinded sample size re-estimation
#> for a stratified design, under equal response probabilities in every stratum
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Planned sizes    : (8, 8) and (8, 8)
#>   Interim sizes    : (4, 4) and (4, 4)
#>   Assumed effect   : 0.3 (RD)
#>   Nominal level    : 0.025
#>   Interim outcomes : 81 pairs (s1, s2)
#>   Above the level  : 10605 outcomes for CRP, 17037 for CRP.total
#> 
#> Largest conditional rejection probability
#>   Quantity   value s1 s2 w1 w2
#>        CRP 0.06250  1  1  1  1
#>  CRP.total 0.05853  2  2  0  0
#> 
#> Type I error rate over 4 points: largest 0.02467 at theta = (0.5, 0.3)
attr(crp, 'TIE')
#>   theta.S1 theta.S2        TIE
#> 1      0.3      0.3 0.02437500
#> 2      0.5      0.3 0.02467373
#> 3      0.3      0.5 0.02467373
#> 4      0.5      0.5 0.02463391
```
