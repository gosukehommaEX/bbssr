# Type I Error Rate of Blinded Sample Size Re-Estimation for a Stratified Design

Calculates the type I error rate of a two-arm trial with a binary
endpoint, stratified into two or more strata, with blinded sample size
re-estimation, over the common response probabilities of the strata
under the null hypothesis, together with that of the fixed design with
the planned group sizes. With two strata the rate is computed exactly
and its largest value over a rectangle of probabilities is certified;
with three or more strata the rate of the re-estimation design is
estimated by simulation.

## Usage

``` r
StratTypeIErrorBSSR(
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
  nsim = 1e+05,
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

- theta:

  Common response probabilities of the strata under the null hypothesis.
  A vector gives, with two strata, the grid of all pairs of its values,
  and with three or more strata the points at which every stratum has
  the same value. A matrix with one column per stratum gives the points
  directly. Default is `seq(0, 1, by = 0.05)`

- maximize:

  How the largest type I error rate is located. With `'certified'` (the
  default with two strata, and available only with two) it is certified
  over the rectangle from the smallest to the largest value of `theta`
  in each stratum, see Details. With `'grid'` (the default with three or
  more strata) it is the largest value at the points of `theta`

- nsim:

  Number of simulated trials at each point with three or more strata.
  Default is `1e5`. Ignored with two strata

- alloc:

  Allocation of the patients to the strata and the groups. `'fixed'`
  (default) uses the group sizes `N1` and `N2`. The others need two
  strata and make the sizes random, see Details: `'stratum'` (stratum
  sizes binomial, permuted blocks within the strata), `'group'` (stratum
  sizes fixed, group sizes binomial) and `'both'` (stratum and group
  sizes binomial)

- block:

  Size of the permuted blocks of `alloc = 'stratum'`. Each block holds
  `block r / (1 + r)` patients of group 1, which must be a whole number.
  Default is 4

- ss.frac:

  Fractions of the strata used by the re-estimation with a random
  allocation: `'observed'` (default), those of the interim patients, or
  `'planned'`, those given by `frac`

- nseq:

  Number of allocation sequences drawn with a random allocation. Default
  is 100

## Value

An object of class `bbssr_stie`, a data frame with one row per point,
containing the probabilities `theta.S1, theta.S2, ...`, the type I error
rates `TIE.BSSR` and `TIE.TRAD` of the designs with re-estimation and
without it, and, with three or more strata or a random allocation, the
Monte Carlo standard error `se.BSSR`. The attribute `max` is a data
frame with the location, the value and, with `maximize = 'certified'`,
the upper bound `bound` of the largest rate of each design, and with a
random allocation the standard error `se` of the rate of the BSSR design
at the point of `theta` with the largest rate.

## Details

Under the null hypothesis both groups of stratum `j` respond with the
probability `theta_j`. Given the interim and the second-stage numbers of
responders of a stratum, the numbers of responders of group 1 are
hypergeometric whatever `theta_j`, so with two strata the type I error
rate is a polynomial in `(theta_1, theta_2)`. Its coefficients in the
products of the Bernstein bases of the two strata are computed exactly,
and the rate is evaluated at the points of `theta` from them. Over a
rectangle the polynomial lies between the smallest and the largest of
its coefficients there, so with `maximize = 'certified'` the rectangle
is halved along its longer side until the largest coefficient of every
piece exceeds the largest value found by at most 1e-12. The largest
value is reported with this upper bound, which the type I error rate
cannot exceed anywhere in the rectangle. The fixed design is certified
in the same way.

With three or more strata the type I error rate of the re-estimation
design is estimated at each point from `nsim` simulated trials, with the
random numbers of R, and that of the fixed design is computed exactly.

With a random allocation (two strata only) the allocation sequences are
drawn as in
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md).
Given a sequence the type I error rate is a polynomial as above, and the
rate of the re-estimation design is the average of these polynomials
over `nseq` sequences, itself a polynomial, whose largest value is
certified as above. The certificate holds for this average, and
`se.BSSR` gives the Monte Carlo standard error of the average at each
point, from the variation between the sequences; the attribute `max`
gives it at the point of `theta` with the largest rate. The rate of the
fixed design with the planned total and the same allocation is exact:
the sum over the stratum and group sizes of their probability times the
polynomial of the fixed design they give.

## See also

[`StratAlphaAdjBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratAlphaAdjBSSR.md),
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
tie <- StratTypeIErrorBSSR(Delta.A = 0.2, N1 = c(20, 20), N2 = c(20, 20), omega = 0.5,
                           r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                           ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.2))
attr(tie, 'max')
#>   Design theta.S1 theta.S2        TIE      bound
#> 1   BSSR      0.5      0.5 0.02508517 0.02508517
#> 2   TRAD      0.5      0.5 0.02833680 0.02833680
```
