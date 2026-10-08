# Adjusted Significance Level of a Blinded Sample Size Re-estimation Design

Lowers the nominal significance level of a two-arm trial with a binary
endpoint and blinded sample size re-estimation (BSSR) until the largest
type I error rate over the common response probability, or over the
pooled response probability on the boundary of a non-inferiority
hypothesis, does not exceed the target level `alpha`, and returns the
adjusted level together with that of the fixed-sample design.

## Usage

``` r
BinaryAlphaAdjBSSR(
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
  n.grid = 100,
  bb.gamma = 0,
  effect = c("RD", "RR", "OR"),
  ss.method = c("exact", "standard", "null.variance", "alternative.variance"),
  ss.Test = Test,
  ss.alpha = alpha,
  rounding = c("group", "friede-kieser", "total", "nearest"),
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50),
  N.min = NULL,
  N.max = NULL,
  n.interim = NULL,
  theta = seq(0, 1, by = 0.005),
  maximize = c("certified", "refined", "grid"),
  adjust = c("test", "both"),
  tol = 1e-08,
  step = 1e-05,
  margin = 0,
  ref.pvalue = FALSE
)
```

## Arguments

- Delta.A:

  Assumed treatment effect, on the scale given by `effect`, used to
  split the blinded pooled proportion into group-specific proportions

- N1:

  Initial sample size of group 1

- N2:

  Initial sample size of group 2

- omega:

  Fraction of the initial sample size observed at the interim analysis.
  The interim size of group 2 is `ceiling(omega N2)` and that of group 1
  is `ceiling(r ceiling(omega N2))`, so the interim analysis keeps the
  allocation ratio. Either `omega` or `n.interim` must be supplied

- r:

  Allocation ratio to group 1

- alpha:

  Target level of significance, which the type I error rate must not
  exceed

- tar.power:

  Target power

- Test:

  Type of statistical test of the final analysis. Options: `'Chisq'`,
  `'Fisher'`, `'Fisher-midP'`, `'Z-pool'`, `'Boschloo'`, `'Blackwelder'`
  or `'Farrington-Manning'`

- restricted:

  Logical. If `TRUE`, the re-estimated sample size is not allowed to
  fall below the initial sample size. Default is `FALSE`

- alternative:

  Direction of the alternative hypothesis. Options: `'greater'`
  (default), `'less'` or `'two.sided'`

- tsmethod:

  Convention used to construct the two-sided version of the conditional
  tests, see
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md).
  Options: `'minlike'` (default), `'central'` or `'blaker'`

- n.grid:

  Number of grid points used to search over the nuisance parameter of
  the unconditional tests. Default is 100

- bb.gamma:

  Confidence level parameter of the Berger-Boos procedure. The default
  of 0 disables the procedure

- effect:

  Scale of `Delta.A` and `Delta.T`. Options: `'RD'` (default) for the
  risk difference `p1 - p2`, `'RR'` for the risk ratio `p1 / p2` or
  `'OR'` for the odds ratio

- ss.method:

  How the sample size is re-estimated from the recovered proportions.
  Options: `'exact'` (default), `'standard'`, `'null.variance'` or
  `'alternative.variance'`, as in
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)

- ss.Test:

  Test whose exact power is used for the re-estimation when
  `ss.method = 'exact'`. Default is `Test`

- ss.alpha:

  Level of significance used for the re-estimation when
  `adjust = 'test'`. Default is `alpha`. It is ignored when
  `adjust = 'both'`

- rounding:

  How the re-estimated sample size is turned into whole numbers.
  Options: `'group'` (default), `'friede-kieser'`, `'total'` or
  `'nearest'`. Only `'group'` is available with `ss.method = 'exact'`.
  See Details

- search:

  How the exact re-estimation chooses the sample size when
  `ss.method = 'exact'`: `'crossing'` (default), `'smallest'` or
  `'stable'`, as in
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).
  Ignored by the other methods

- search.limit:

  Limit of the size of group 2 examined by `search = 'stable'`, as in
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).
  Default is `c(2, 50)`

- N.min:

  Lower bound on the final total sample size, or `NULL` (default) for
  none beyond the interim total. It can be used to keep the patients who
  are already enrolled but not yet evaluated at the interim analysis

- N.max:

  Upper bound on the final total sample size, or `NULL` (default) for
  none

- n.interim:

  Interim sample sizes of group 1 and group 2, as a vector of length
  two. An alternative to `omega`

- theta:

  Grid of common response probabilities at which the type I error rate
  is evaluated. Default is `seq(0, 1, by = 0.005)`. With a non-zero
  `margin`, `theta` is the pooled response probability
  `(r p1 + p2) / (1 + r)` on the boundary of the null hypothesis, see
  Details

- maximize:

  How the largest type I error rate at a level is located, as in
  [`BinaryTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md).
  With `'certified'` (default) the adjusted level controls the type I
  error rate over the whole interval from the smallest to the largest
  value of `theta`, restricted with a margin to the values at which both
  response probabilities on the null boundary lie in the unit interval,
  see Details. With `'refined'` and `'grid'` the control is assessed on
  the grid

- adjust:

  Which part of the design uses the adjusted level. `'test'` (default)
  applies it to the final analysis only, so the re-estimation keeps
  using `ss.alpha`. `'both'` applies it to the re-estimation as well

- tol:

  Tolerance of the bisection, relative to `alpha`, a single value in (0,
  1). The bisection is used when `adjust = 'test'` and, whatever
  `adjust`, for the fixed-sample design. Default is `1e-8`

- step:

  Decrement of the level in the search used when `adjust = 'both'`.
  Default is `1e-5`

- margin:

  Non-inferiority margin on the scale of the risk difference. The
  default of 0 gives a test of superiority. A value other than 0 tests
  the null hypothesis `p1 - p2 <= -margin` against `p1 - p2 > -margin`
  when `alternative` is `'greater'`, and `p1 - p2 >= margin` against
  `p1 - p2 < margin` when it is `'less'`. It requires
  `Test = 'Blackwelder'` or `'Farrington-Manning'`, see
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md).
  A negative value tests for superiority by more than its absolute
  value. With a value other than 0 the assumed and the true effects are
  risk differences (`effect = 'RD'`)

- ref.pvalue:

  Logical. If `TRUE`, the maximization over the nuisance parameter of
  the unconditional tests is refined between the grid points, see
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md).
  Default is `FALSE`. It applies to the final analysis, to the
  fixed-sample comparator and to the exact re-estimation

## Value

An object of class `bbssr_alphaadj`, a data frame with one row for the
BSSR design and one for the fixed-sample design containing:

- Design:

  `'BSSR'` or `'TRAD'`

- alpha:

  Target level

- max.TIE:

  Largest type I error rate at the nominal level `alpha`

- alpha.adj:

  Adjusted nominal level, see Details. It is `alpha` when the nominal
  level already controls the type I error rate, and 0, with a warning,
  when no level examined does

- max.TIE.adj:

  Largest type I error rate at the adjusted level, 0 when `alpha.adj` is
  0

- theta.adj:

  Common response probability, or pooled response probability on the
  null boundary, at which `max.TIE.adj` occurs, `NA` when `alpha.adj` is
  0

- max.TIE.bound:

  Upper bound of the type I error rate at the nominal level with
  `maximize = 'certified'`, otherwise `NA`

- max.TIE.adj.bound:

  Upper bound of the type I error rate at the adjusted level with
  `maximize = 'certified'`, otherwise `NA`

## Details

Following Kieser and Friede (2000) and Friede and Kieser (2004), the
nominal level of the final analysis is lowered until the largest type I
error rate over the common response probability, evaluated as in
[`BinaryTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md),
does not exceed `alpha`.

With `adjust = 'test'` the final sample sizes do not depend on the
nominal level, and the rejection regions shrink as the level decreases,
so the largest type I error rate is a non-decreasing step function of
the level. The bisection locates, to within `tol * alpha`, the level
above which the largest rate exceeds `alpha`. The package rejects a
p-value at a level when it is below the level by more than a tolerance
of about 1.49e-8, so every level between the largest p-value rejected at
the level found plus this tolerance and the smallest p-value not
rejected gives the same rejection regions, also when a p-value is
rejected if it is below the level or at most equal to it. The reported
`alpha.adj` is the largest value with at most six significant digits in
this range, or with more digits if there is none, so the printed value
can be used as it is. The same applies to the fixed-sample design.

With `adjust = 'both'` the re-estimated sample sizes change with the
level as well, so the type I error rate need not be monotone in the
level. The level is then lowered from `alpha` in steps of `step` until
the type I error rate is controlled. The result is the first level among
`alpha - step`, `alpha - 2 * step`, and so on, at which the type I error
rate is controlled, and it can change with `step`. The levels that
control the type I error rate need not form an interval, so a level
below the result can fail, and a level above it that is not on this grid
can pass. This search re-estimates the sample size at every step and can
take much longer than the bisection. The fixed-sample design always uses
the bisection, since its sample size does not depend on the level. The
script `reproduce-published.R` in the folder
`system.file('reproduce', package = 'bbssr')` reproduces the published
values of Section 5 of Friede and Kieser (2004) and of Example 21.1 of
Kieser (2020) only when the re-estimation also uses the adjusted level,
as with `adjust = 'both'`.

The type I error rate can also be controlled within a confidence
interval for the common response probability, or for the pooled response
probability on the null boundary, instead of over its whole range
(Kieser, 2020, Section 22.1.3). A value `gamma` between 0 and `alpha` is
fixed in advance, and the final analysis uses a level at which the
largest type I error rate over the `1 - gamma` confidence interval
computed from the data of the completed trial does not exceed
`alpha - gamma`. The type I error rate is then controlled at `alpha`.
This level is found when `theta` spans the interval and `alpha` is set
to `alpha - gamma`. The interval is known only at the end of the trial,
so the re-estimation keeps the nominal level, with `adjust = 'test'` and
`ss.alpha` set to the nominal level. Neither search examines levels
above `alpha - gamma`. The script `reproduce-published.R` reproduces the
adjusted level of Example 23.1 of Kieser (2020), which uses the
Clopper-Pearson interval with `gamma = 1e-4`.

With `maximize = 'certified'` the levels are assessed on the grid with
the refinement of `maximize = 'refined'`, and the level found is then
certified as in
[`BinaryTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md):
it is accepted only if the upper bound of the type I error rate over the
interval from the smallest to the largest value of `theta` does not
exceed `alpha`. If the bound exceeds `alpha`, the grid has missed a
value of `theta` at which the level fails, or the largest rate lies
within 1e-12 below `alpha`, and the search continues below that level.
The bisection then certifies every level that passes the assessment on
the grid, and the search of `adjust = 'both'` certifies every such level
from the start. The adjusted level therefore controls the type I error
rate over the whole interval. If a certification stops at its limit of
1e5 subdivisions, with a warning, the bound can exceed the largest rate
by more than 1e-12, and the level found can be lower than necessary.

In both searches the type I error rate at a new level is first evaluated
at the single value of `theta` where the largest rate of the last level
that failed was found. If it exceeds `alpha` there, the level fails
without the evaluation over the whole grid. The largest rate over
`theta` is never below the rate at any one value, so such a level also
fails on the largest rate. With `maximize = 'grid'` this value is a grid
point, and every decision is that of the evaluation on the grid. With
`'refined'` and `'certified'` it can lie between the grid points, where
the grid can miss a rate above `alpha`, so a level can fail that the
evaluation on the grid would accept. This can only lower the level
found, and with `'certified'` the level found is certified in either
case.

## References

Kieser M, Friede T (2000). Re-calculating the sample size in internal
pilot study designs with control of the type I error rate. *Statistics
in Medicine*, 19(7), 901-911.

Friede T, Kieser M (2004). Sample size recalculation for binary data in
internal pilot study designs. *Pharmaceutical Statistics*, 3(4),
269-279.

Kieser M (2020). *Methods and Applications of Sample Size Calculation
and Recalculation in Clinical Trials*. Springer, Cham.

## See also

[`BinaryTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
# \donttest{
BinaryAlphaAdjBSSR(
  Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
  theta = seq(0.01, 0.99, by = 0.01)
)
#> Adjusted significance level for blinded sample size re-estimation
#> 
#>   Test            : Chisq
#>   Alternative     : greater
#>   Adjusted part   : final analysis only
#>   Maximum         : certified over theta in [0.01, 0.99]
#>   Target level    : 0.025
#> 
#>        Design   max.TIE alpha.adj max.TIE.adj
#>          BSSR 0.0272864 0.0204375   0.0237418
#>  Fixed sample 0.0293761 0.0207700   0.0247660
# }
```
