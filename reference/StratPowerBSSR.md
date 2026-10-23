# Power of Blinded Sample Size Re-Estimation for a Stratified Design

Calculates the power and the distribution of the final total sample size
of a two-arm trial with a binary endpoint, stratified into two or more
strata, with blinded sample size re-estimation, together with the power
of the fixed design with the planned group sizes. The sample size is
re-estimated at an interim analysis from the pooled numbers of
responders of the strata, without the treatment groups, and the final
analysis is a test of
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md)
with all patients. With two strata the power is computed exactly, and
with three or more strata by simulation.

## Usage

``` r
StratPowerBSSR(
  p,
  Delta.A,
  Delta.T,
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
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50),
  N.min = NULL,
  N.max = NULL,
  n.interim = NULL,
  frac = NULL,
  nsim = 1e+05
)
```

## Arguments

- p:

  True pooled response probabilities of the strata: a vector with one
  element per stratum, or a matrix with one row per scenario and one
  column per stratum

- Delta.A:

  Assumed treatment effect, on the scale given by `effect`, common to
  the strata and used to divide the blinded estimate of the pooled
  response probability of each stratum between the groups at the
  re-estimation

- Delta.T:

  True treatment effect, on the scale given by `effect`: a single value
  common to the strata, or one value per stratum

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
  one-sided levels up to 0.078) or `'Strat-Fisher'` (the stratified
  exact conditional test of Jung (2014)). See Details

- restricted:

  Logical. If `TRUE`, the final group sizes are at least the planned
  sizes. Default is `FALSE`

- alternative:

  Direction of the alternative hypothesis. Options: `'greater'`
  (default), under which group 1 responds more often than group 2,
  `'less'` or `'two.sided'`

- tsmethod:

  Convention of the two-sided p-value of `'Strat-Fisher'`, as for the
  Fisher exact test of
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

- n.interim:

  Interim group sizes, as a matrix with two rows (group 1 and group 2)
  and one column per stratum. Supply either `omega` or `n.interim`

- frac:

  Fractions of the strata used to divide the re-estimated total sample
  size between the strata, see
  [`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md).
  Default is `NULL`, which takes the fractions of the planned design,
  `(N1 + N2) / sum(N1 + N2)`

- nsim:

  Number of simulated trials for each scenario with three or more
  strata. Default is `1e5`. Ignored with two strata

## Value

An object of class `bbssr_spowerbssr`, a data frame with one row for
each scenario in which both response probabilities of every stratum lie
in the unit interval, containing:

- p.S1, p.S2, ...:

  True pooled response probabilities of the strata

- p1.S1, ..., p2.S1, ...:

  True response probabilities of the two groups

- power.BSSR:

  Power of the design with re-estimation

- power.TRAD:

  Power of the fixed design with the planned group sizes

- E.N:

  Expected final total sample size

- se.BSSR:

  Monte Carlo standard error of `power.BSSR`, with three or more strata
  only

The attribute `reestimation` holds, with two strata, the final group
sizes for every pair of interim pooled numbers of responders, `N.dist`
the distribution of the final total sample size in each scenario, and
the other attributes the design settings.

## Details

In stratum `j` the blinded estimate of the pooled response probability
is the number `s_j` of interim responders divided by the interim size of
the stratum. It is divided between the groups with the assumed effect
`Delta.A` and the allocation ratio of the stratum, as in
[`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md),
and the total sample size is re-estimated from these probabilities by
the rule `ss.method`. The total is divided between the strata and the
groups as in
[`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md),
and the final group sizes are at least the interim sizes, at least the
planned sizes when `restricted = TRUE`, with a total of at least `N.min`
and at most `N.max`. When the probabilities of every stratum show no
effect in the direction of the alternative, which happens with
`effect = 'RR'` when no interim patient responds and with
`effect = 'OR'` when in each stratum no interim patient or every interim
patient responds, the planned group sizes are kept. A stratum whose
pooled estimate is close to 0 with `effect = 'RR'`, or close to 0 or 1
with `effect = 'OR'`, has planning probabilities close to each other, so
the re-estimated size can be very large, and the time of the exact
computation grows with it; an upper bound `N.max` keeps both in check.

With two strata the final group sizes depend on the interim data only
through `(s_1, s_2)`. Given `s_j`, the number of interim responders of
group 1 in stratum `j` follows from the two binomial distributions of
the stratum, and the second-stage numbers of responders are binomial
given the final sizes. The power is the sum over `(s_1, s_2)` of the
probability of the interim outcome times the conditional probability of
rejection, computed exactly by enumeration in compiled code.

With three or more strata the number of interim outcomes grows as a
product over the strata, and the power is estimated from `nsim`
simulated trials, with the random numbers of R, so that
[`set.seed()`](https://rdrr.io/r/base/Random.html) makes the results
reproducible. The power of the fixed design is exact for any number of
strata.

## References

Woolson RF, Bean JA, Rojas PB (1986). Sample size for case-control
studies using Cochran's statistic. *Biometrics*, 42(4), 927-932.

Nam J (1992). Sample size determination for case-control studies and the
comparison of stratified and unstratified analyses. *Biometrics*, 48(2),
389-395.

Jung SH (2014). Stratified Fisher's exact test and its sample size
calculation. *Biometrical Journal*, 56(1), 129-140.

## See also

[`StratTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md),
[`StratBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
res <- StratPowerBSSR(p = rbind(c(0.3, 0.5), c(0.4, 0.6)), Delta.A = 0.2,
                      Delta.T = 0.2, N1 = c(30, 30), N2 = c(30, 30), omega = 0.5,
                      r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                      ss.method = 'standard')
print(res)
#> Blinded sample size re-estimation for a stratified design
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Re-estimation    : normal approximation (unrestricted)
#>   Planned sizes    : (30, 30) and (30, 30)
#>   Interim sizes    : (15, 15) and (15, 15)
#>   Effect           : RD, assumed 0.2, true 0.2
#>   Alpha            : 0.025, target power 0.8
#> 
#>  p.S1 p.S2 p1.S1 p1.S2 p2.S1 p2.S2 power.BSSR power.TRAD    E.N
#>   0.3  0.5   0.4   0.6   0.2   0.4     0.7892     0.6274 174.55
#>   0.4  0.6   0.5   0.7   0.3   0.5     0.7888     0.6162 181.99
```
