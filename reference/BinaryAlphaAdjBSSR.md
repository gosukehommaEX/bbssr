# Adjusted Significance Level of a Blinded Sample Size Re-estimation Design

Finds the largest nominal significance level at which a two-arm trial
with a binary endpoint and blinded sample size re-estimation (BSSR)
keeps its type I error rate at or below the target level `alpha` for
every common response probability, together with the corresponding level
of the fixed-sample design.

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
  tsmethod = c("minlike", "central"),
  n.grid = 100,
  bb.gamma = 0,
  effect = c("RD", "RR", "OR"),
  ss.method = c("exact", "standard", "null.variance", "alternative.variance"),
  ss.Test = Test,
  ss.alpha = alpha,
  rounding = c("group", "friede-kieser", "total", "nearest"),
  N.min = NULL,
  N.max = NULL,
  n.interim = NULL,
  theta = seq(0.005, 0.995, by = 0.005),
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
  tests. Options: `'minlike'` (default) or `'central'`

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
  is evaluated. Default is `seq(0.005, 0.995, by = 0.005)`. With a
  non-zero `margin`, `theta` is the pooled response probability
  `(r p1 + p2) / (1 + r)` on the boundary of the null hypothesis, see
  Details

- adjust:

  Which part of the design uses the adjusted level. `'test'` (default)
  applies it to the final analysis only, so the re-estimation keeps
  using `ss.alpha`. `'both'` applies it to the re-estimation as well

- tol:

  Tolerance of the bisection used when `adjust = 'test'`, relative to
  `alpha`. Default is `1e-8`

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

  Adjusted nominal level

- max.TIE.adj:

  Largest type I error rate at the adjusted level

- theta.adj:

  Common response probability, or pooled response probability on the
  null boundary, at which `max.TIE.adj` occurs

## Details

Following Kieser and Friede (2000) and Friede and Kieser (2004), the
nominal level of the final analysis is lowered until the largest type I
error rate over the common response probability, evaluated as in
[`BinaryTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md),
does not exceed `alpha`.

With `adjust = 'test'` the final sample sizes do not depend on the
nominal level, and the rejection regions shrink as the level decreases,
so the largest type I error rate is a non-decreasing step function of
the level. The adjusted level is then found by bisection to within
`tol * alpha`, and the reported value is the largest level examined that
controls the type I error rate.

With `adjust = 'both'` the re-estimated sample sizes change with the
level as well, so the type I error rate need not be monotone in the
level. The level is then lowered from `alpha` in steps of `step` until
the type I error rate is controlled. This search re-estimates the sample
size at every step and can take much longer than the bisection. The
fixed-sample design always uses the bisection, since its sample size
does not depend on the level.

In both searches the type I error rate at a new level is first evaluated
at the single grid point where the largest rate of the last level that
failed was found. If it exceeds `alpha` there, the level fails without
the evaluation over the whole grid and the refinement. The largest rate
over `theta` is never below the rate at a grid point, so every decision,
and hence the result, is that of the full evaluation.

## References

Kieser M, Friede T (2000). Re-calculating the sample size in internal
pilot study designs with control of the type I error rate. *Statistics
in Medicine*, 19(7), 901-911.

Friede T, Kieser M (2004). Sample size recalculation for binary data in
internal pilot study designs. *Pharmaceutical Statistics*, 3(4),
269-279.

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
#>   Target level    : 0.025
#> 
#>        Design   max.TIE alpha.adj max.TIE.adj
#>          BSSR 0.0272864 0.0204376   0.0237418
#>  Fixed sample 0.0293761 0.0207700   0.0247660
# }
```
