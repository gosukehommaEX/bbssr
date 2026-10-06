# Sample Size Re-estimation from Observed Blinded Interim Data

Re-estimates the sample size of an ongoing two-arm trial with a binary
endpoint from the blinded data available at an interim analysis, and
reports how many patients still have to be enrolled in each group during
the second stage. Only the total number of patients and the total number
of responders are required, so the treatment allocation remains
concealed.

## Usage

``` r
BinaryBSSR(
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
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker"),
  n.grid = 100,
  bb.gamma = 0,
  effect = c("RD", "RR", "OR"),
  ss.method = c("exact", "standard", "null.variance", "alternative.variance"),
  ss.Test = Test,
  ss.alpha = alpha,
  rounding = c("group", "friede-kieser", "total", "nearest"),
  N.min = NULL,
  N.max = NULL,
  margin = 0,
  ref.pvalue = FALSE
)
```

## Arguments

- n1:

  Number of patients of group 1 observed at the interim analysis

- n2:

  Number of patients of group 2 observed at the interim analysis

- S:

  Total number of responders observed at the interim analysis, pooled
  over both groups

- Delta.A:

  Assumed treatment effect, on the scale given by `effect`, used to
  split the blinded pooled proportion into group-specific proportions

- r:

  Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2
  = r:1, r \> 0)

- alpha:

  Level of significance of the final analysis, for the alternative
  specified by `alternative`

- tar.power:

  Target power

- Test:

  Type of statistical test of the final analysis. Options: `'Chisq'`,
  `'Fisher'`, `'Fisher-midP'`, `'Z-pool'`, `'Boschloo'`, `'Blackwelder'`
  or `'Farrington-Manning'`

- restricted:

  Logical. If `TRUE`, the re-estimated sample size is not allowed to
  fall below the planned sample size given by `N1` and `N2`. Default is
  `FALSE`

- N1:

  Planned sample size of group 1. Required when `restricted` is `TRUE`,
  and when the recovered proportions coincide

- N2:

  Planned sample size of group 2. Required when `restricted` is `TRUE`,
  and when the recovered proportions coincide

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

  Scale of `Delta.A`. Options: `'RD'` (default), `'RR'` or `'OR'`, as in
  [`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)

- ss.method:

  How the sample size is re-estimated. Options: `'exact'` (default),
  `'standard'`, `'null.variance'` or `'alternative.variance'`, as in
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)

- ss.Test:

  Test whose exact power is used for the re-estimation when
  `ss.method = 'exact'`. Default is `Test`

- ss.alpha:

  Level of significance used for the re-estimation. Default is `alpha`

- rounding:

  How the re-estimated sample size is turned into whole numbers.
  Options: `'group'` (default), `'friede-kieser'`, `'total'` or
  `'nearest'`, as in
  [`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)

- N.min:

  Lower bound on the final total sample size, or `NULL` (default). It
  can be used to keep the patients who are already enrolled but not yet
  evaluated

- N.max:

  Upper bound on the final total sample size, or `NULL` (default)

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
  Default is `FALSE`

## Value

An object of class `bbssr_bssr`, a data frame with one row containing:

- n1:

  Interim sample size of group 1

- n2:

  Interim sample size of group 2

- n:

  Total interim sample size

- S:

  Total number of interim responders

- hat.p:

  Blinded estimate of the pooled response probability

- hat.p1:

  Recovered response probability of group 1

- hat.p2:

  Recovered response probability of group 2

- N1.re:

  Re-estimated total sample size of group 1

- N2.re:

  Re-estimated total sample size of group 2

- N.re:

  Re-estimated total sample size

- n1.stage2:

  Number of additional patients to enrol in group 1

- n2.stage2:

  Number of additional patients to enrol in group 2

- n.stage2:

  Total number of additional patients to enrol

- N1.final:

  Final sample size of group 1

- N2.final:

  Final sample size of group 2

- N.final:

  Final total sample size

- Power:

  Exact power at the final sample size under the recovered proportions

## Details

The blinded estimate of the pooled response probability is
`hat.p = S / (n1 + n2)`. For the risk difference, group-specific
proportions are recovered as `hat.p1 = hat.p + Delta.A / (1 + r)` and
`hat.p2 = hat.p - r Delta.A / (1 + r)`, truncated to the unit interval.
For the risk ratio and the odds ratio they are the proportions with the
pooled value `hat.p` and the ratio `Delta.A`. The sample size is then
re-estimated as in
[`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).

Under the unrestricted rule and `rounding = 'group'`, the final size of
group 2 is the larger of the re-estimated size and what has already been
observed. Under the restricted rule it is raised to the planned size
first, so the trial can only grow, and `N.min` and `N.max` add a further
lower and upper bound on the total. The final size of group 1 is then
`ceiling(r N2.final)`, so an imbalance already present at the interim is
corrected by the remaining enrolment instead of being carried forward.
Neither second-stage size is ever negative. The other rounding rules are
described in
[`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md).

While
[`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
evaluates the operating characteristics of a BSSR design at the planning
stage, this function is applied once, to the data of a trial that is
under way.

## See also

[`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md),
[`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
# Interim data: 20 patients per group, 11 responders in total
BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
           alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
#> Blinded sample size re-estimation from interim data
#> 
#>   Test             : Chisq
#>   Alternative      : greater
#>   Design rule      : unrestricted
#>   Assumed effect   : 0.3
#>   Alpha            : 0.025, target power 0.8
#> 
#> Interim data
#>   Patients         : n1 = 20, n2 = 20, total n = 40
#>   Responders       : S = 11, blinded pooled rate = 0.275
#>   Recovered rates  : hat.p1 = 0.425, hat.p2 = 0.125
#> 
#> Re-estimation
#>   Required total   : N1 = 33, N2 = 33, total N = 66
#>   Still to enrol   : group 1 = 13, group 2 = 13, total = 26
#>   Final size       : N1 = 33, N2 = 33, total N = 66
#>   Power at final N : 0.8136

# \donttest{
# Restricted rule with a planned sample size of 40 per group
BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
           alpha = 0.025, tar.power = 0.8, Test = 'Boschloo',
           restricted = TRUE, N1 = 40, N2 = 40)
#> Blinded sample size re-estimation from interim data
#> 
#>   Test             : Boschloo
#>   Alternative      : greater
#>   Design rule      : restricted
#>   Assumed effect   : 0.3
#>   Alpha            : 0.025, target power 0.8
#> 
#> Interim data
#>   Patients         : n1 = 20, n2 = 20, total n = 40
#>   Responders       : S = 11, blinded pooled rate = 0.275
#>   Recovered rates  : hat.p1 = 0.425, hat.p2 = 0.125
#> 
#> Re-estimation
#>   Required total   : N1 = 34, N2 = 34, total N = 68
#>   Still to enrol   : group 1 = 20, group 2 = 20, total = 40
#>   Final size       : N1 = 40, N2 = 40, total N = 80
#>   Power at final N : 0.8602
# }
```
