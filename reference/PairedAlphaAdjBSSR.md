# Adjusted Significance Level of Blinded Sample Size Re-Estimation for Paired Binary Data

Lowers the nominal significance level of a matched-pairs design with
blinded sample size re-estimation until the largest type I error rate
over the probability of a discordant pair on the boundary of the null
hypothesis does not exceed the target level `alpha`, and returns the
adjusted level together with that of the fixed design.

## Usage

``` r
PairedAlphaAdjBSSR(
  Delta.A,
  N,
  omega = NULL,
  alpha,
  tar.power,
  Test,
  restricted = FALSE,
  alternative = c("greater", "less", "two.sided", "equivalence"),
  n.grid = 100,
  bb.gamma = 0,
  ss.method = c("exact", "standard", "conditional.power"),
  ss.Test = Test,
  ss.alpha = alpha,
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50),
  N.min = NULL,
  N.max = NULL,
  n.interim = NULL,
  psi = seq(0, 1, by = 0.01),
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

  Assumed difference `p10 - p01` used in the re-estimation

- N:

  Planned number of pairs, which is the number of pairs of the fixed
  design

- omega:

  Fraction of `N` at the interim analysis, rounded up. Supply either
  `omega` or `n.interim`

- alpha:

  Target level of significance, which the type I error rate must not
  exceed

- tar.power:

  Target power of the re-estimation

- Test:

  Type of statistical test. Options: `'McNemar'`, `'McNemar-cc'`,
  `'McNemar-exact'`, `'McNemar-midP'`, `'McNemar-uncond'` or
  `'Paired-Score'`

- restricted:

  Logical. If `TRUE`, the final number of pairs is at least `N`. Default
  is `FALSE`

- alternative:

  Direction of the alternative hypothesis. Options: `'greater'`
  (default) for `p10 - p01 > -margin`, `'less'` for
  `p10 - p01 < margin`, `'two.sided'` for `p10 != p01`, or
  `'equivalence'` for `-lower < p10 - p01 < upper`

- n.grid:

  Number of grid points used to search over the probability of a
  discordant pair in the test `'McNemar-uncond'`. Default is 100

- bb.gamma:

  Confidence level parameter of the Berger-Boos procedure for the test
  `'McNemar-uncond'`. The default of 0 disables the procedure

- ss.method:

  Re-estimation rule. `'exact'` (default) and `'standard'` re-estimate
  the number of pairs of
  [`PairedSampleSize`](https://gosukehommaex.github.io/bbssr/reference/PairedSampleSize.md)
  with the method of that name. `'conditional.power'` takes the smallest
  number of pairs at which the conditional power given the interim
  number of discordant pairs attains `tar.power`, the rule of Schepers
  et al. (2026); it requires `N.max`. See Details

- ss.Test:

  Test used in the re-estimation. Default is `Test`

- ss.alpha:

  Level of significance used for the re-estimation when
  `adjust = 'test'`. Default is `alpha`. It is ignored when
  `adjust = 'both'`

- search:

  Search of the exact re-estimation, see
  [`PairedSampleSize`](https://gosukehommaex.github.io/bbssr/reference/PairedSampleSize.md);
  it does not go beyond `N.max`

- search.limit:

  Limit of the search `'stable'`, see
  [`PairedSampleSize`](https://gosukehommaex.github.io/bbssr/reference/PairedSampleSize.md)

- N.min:

  Smallest final number of pairs, or `NULL`

- N.max:

  Largest final number of pairs, or `NULL`

- n.interim:

  Number of pairs at the interim analysis. Supply either `omega` or
  `n.interim`

- psi:

  Values of the probability `p10 + p01` of a discordant pair at which
  the rate is evaluated. Values that do not lie on the boundary of the
  null hypothesis, below `|margin|`, are dropped. Default is
  `seq(0, 1, by = 0.01)`

- maximize:

  How the largest type I error rate at a level is located, as in
  [`PairedTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedTypeIErrorBSSR.md).
  With `'certified'` (default) the adjusted level controls the type I
  error rate over the whole interval spanned by the values of `psi` on
  the boundary. With `'refined'` and `'grid'` the control is assessed on
  the grid

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

- margin:

  Margin on the scale of the difference `p10 - p01`. The default of 0
  gives a test of superiority, and a value other than 0 tests the null
  hypothesis `p10 - p01 <= -margin` against `p10 - p01 > -margin` when
  `alternative` is `'greater'`, and `p10 - p01 >= margin` against
  `p10 - p01 < margin` when it is `'less'`. A negative value tests for
  superiority by more than its absolute value. A margin other than 0
  requires `Test = 'Paired-Score'` and a one-sided alternative. For
  `alternative = 'equivalence'` the margin is the pair
  `c(lower, upper)`, both in (0, 1), or a single value used for both,
  and the null hypothesis `p10 - p01 <= -lower` or `p10 - p01 >= upper`
  is tested

- ref.pvalue:

  Logical. If `TRUE`, the maximization over the probability of a
  discordant pair in the test `'McNemar-uncond'` is refined between the
  grid points, see Details. Default is `FALSE`

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

- psi.adj:

  Probability of a discordant pair at which `max.TIE.adj` occurs, `NA`
  when `alpha.adj` is 0

- max.TIE.bound:

  Upper bound of the type I error rate at the nominal level with
  `maximize = 'certified'`, otherwise `NA`

- max.TIE.adj.bound:

  Upper bound of the type I error rate at the adjusted level with
  `maximize = 'certified'`, otherwise `NA`

For `alternative = 'equivalence'` the column `boundary.adj` gives the
boundary, `'lower'` or `'upper'`, on which `max.TIE.adj` occurs.

## Details

The searches are those of
[`BinaryAlphaAdjBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md).
With `adjust = 'test'` the final numbers of pairs do not depend on the
nominal level of the final analysis, the largest type I error rate is a
non-decreasing function of the level, and the level is found by
bisection to within `tol * alpha`; the reported level is the largest
value with at most six significant digits that gives the same rejection
regions. With `adjust = 'both'` the level is lowered from `alpha` in
steps of `step` until the type I error rate is controlled. The fixed
design always uses the bisection. With `maximize = 'certified'` a level
is accepted only when the upper bound of the type I error rate over the
interval, computed as in
[`PairedTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedTypeIErrorBSSR.md),
does not exceed `alpha`.

The exact conditional McNemar test keeps its level whatever blinded rule
is used (Schepers et al., 2026), so it needs no adjustment.

## References

Kieser M, Friede T (2000). Re-calculating the sample size in internal
pilot study designs with control of the type I error rate. *Statistics
in Medicine*, 19(7), 901-911.

Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026).
Blinded sample size review for McNemar's test based on primary and
surrogate endpoints. arXiv:2608.17784.

## Examples

``` r
# \donttest{
PairedAlphaAdjBSSR(Delta.A = 0.2, N = 40, n.interim = 20, alpha = 0.05,
                   tar.power = 0.8, Test = 'McNemar', alternative = 'two.sided',
                   ss.method = 'standard', psi = seq(0.05, 0.95, by = 0.05))
#> Adjusted significance level for blinded sample size re-estimation
#> 
#>   Test            : McNemar
#>   Alternative     : two.sided
#>   Adjusted part   : final analysis only
#>   Maximum         : certified over psi in [0.05, 0.95]
#>   Target level    : 0.05
#> 
#>        Design   max.TIE alpha.adj max.TIE.adj
#>          BSSR 0.0505829 0.0496472   0.0495041
#>  Fixed sample 0.0538321 0.0481934   0.0485610
# }
```
