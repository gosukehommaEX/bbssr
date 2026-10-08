# Sample Size Calculation for Two-Arm Trials with Binary Endpoints

Calculates the required sample size for two-arm trials with binary
endpoints. Seven tests are supported (see
[`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)),
each of which can be applied with a one-sided or a two-sided
alternative. The sample size can be obtained from the exact power of the
selected test or from the normal approximation.

## Usage

``` r
BinarySampleSize(
  p1,
  p2,
  r,
  alpha,
  tar.power,
  Test,
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker"),
  n.grid = 100,
  bb.gamma = 0,
  method = c("exact", "standard", "null.variance", "alternative.variance"),
  rounding = c("group", "friede-kieser", "total", "nearest"),
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50),
  margin = 0,
  ref.pvalue = FALSE,
  margin.scale = c("RD", "RR")
)
```

## Arguments

- p1:

  True probability of responders for group 1

- p2:

  True probability of responders for group 2

- r:

  Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2
  = r:1, r \> 0)

- alpha:

  Level of significance for the alternative specified by `alternative`

- tar.power:

  Target power

- Test:

  Type of statistical test. Options: `'Chisq'`, `'Fisher'`,
  `'Fisher-midP'`, `'Z-pool'`, `'Boschloo'`, `'Blackwelder'` or
  `'Farrington-Manning'`

- alternative:

  Direction of the alternative hypothesis. Options: `'greater'`
  (default), which requires `p1 - p2 > -margin`, `'less'`, which
  requires `p1 - p2 < margin`, or `'two.sided'`, which requires
  `p1 != p2`. With `margin.scale = 'RR'`, `'greater'` requires
  `p1 > margin p2` and `'less'` requires `p1 < margin p2`

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

- method:

  How the sample size is obtained. `'exact'` (default) steps from the
  normal approximation, one patient at a time, to the sample size at
  which the exact power of `Test` first attains the target, see Details.
  `'standard'` uses the normal approximation with the variance under the
  null hypothesis for the significance term and the variance under the
  alternative for the power term, formula (21.3) of Kieser (2020).
  `'null.variance'` uses the variance under the null hypothesis for both
  terms, formula (1) of Friede and Kieser (2004).
  `'alternative.variance'` uses the variance under the alternative for
  both terms, as in Blackwelder (1982). See Details for a non-zero
  `margin` and for a margin on the scale of the risk ratio

- rounding:

  How an unrounded sample size from the normal approximation is turned
  into whole numbers. `'group'` (default) rounds up the size of group 2
  and gives group 1 `ceiling(r N2)` patients. `'friede-kieser'` rounds
  up the two group sizes separately, as in Friede and Kieser (2004).
  `'total'` rounds up the total and gives group 2 `floor(N / (1 + r))`
  patients. `'nearest'` rounds the two group sizes separately to the
  nearest whole number, as in Farrington and Manning (1990). Only
  `'group'` is available with `method = 'exact'`

- search:

  How the search of `method = 'exact'` chooses the size of group 2.
  `'crossing'` (default) steps from the normal approximation to a size
  that attains the target power while the size one unit smaller does
  not. `'smallest'` returns the smallest size that attains the target
  power. `'stable'` returns the smallest size from which every size up
  to a limit attains it. See Details. Ignored by the other methods

- search.limit:

  Two numbers `a` and `b` giving the limit
  `max(ceiling(a n0), ceiling(n0 + b))` of the size of group 2 examined
  by `search = 'stable'`, where `n0` is the size of group 2 from the
  normal approximation, rounded up. `a` must be at least 1 and `b`
  non-negative. Default is `c(2, 50)`

- margin:

  Non-inferiority margin, on the scale given by `margin.scale`. On the
  scale of the risk difference the default of 0 gives a test of
  superiority, and a value other than 0 tests the null hypothesis
  `p1 - p2 <= -margin` against `p1 - p2 > -margin` when `alternative` is
  `'greater'`, and `p1 - p2 >= margin` against `p1 - p2 < margin` when
  it is `'less'`. A negative value tests for superiority by more than
  its absolute value. On the scale of the risk ratio the margin is a
  positive ratio `R0`, and the null hypothesis is `p1 / p2 <= R0`
  against `p1 / p2 > R0` when `alternative` is `'greater'`, and
  `p1 / p2 >= R0` against `p1 / p2 < R0` when it is `'less'`. A margin
  other than 0, and any margin on the scale of the risk ratio, requires
  `Test = 'Blackwelder'` or `'Farrington-Manning'` and a one-sided
  alternative, see
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)

- ref.pvalue:

  Logical. If `TRUE`, the maximization over the nuisance parameter of
  the unconditional tests is refined between the grid points, see
  [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md).
  Default is `FALSE`

- margin.scale:

  Scale of `margin`. Options: `'RD'` (default) for the risk difference
  `p1 - p2` or `'RR'` for the risk ratio `p1 / p2`

## Value

An object of class `bbssr_samplesize`, a data frame with one row
containing:

- p1:

  True probability of responders for group 1

- p2:

  True probability of responders for group 2

- r:

  Allocation ratio to group 1

- alpha:

  Level of significance

- tar.power:

  Target power

- Test:

  Name of the statistical test

- alternative:

  Direction of the alternative hypothesis

- Power:

  Exact power of `Test` at the selected sample size

- N1:

  Required sample size of group 1

- N2:

  Required sample size of group 2

- N:

  Total required sample size

The attribute `search` holds the search, and the attribute
`search.limit` the limit of the size of group 2 examined by
`search = 'stable'` (otherwise `NA`).

## Details

With the default `search = 'crossing'` the calculation uses a three-step
approach:

1.  Calculate an initial sample size from the normal approximation to
    the chi-squared test

2.  Evaluate the exact power at the initial sample size

3.  If the exact power attains the target power, lower the size of group
    2 one unit at a time as long as the power still attains it;
    otherwise raise it one unit at a time until the power attains it

The normal approximation of the first step uses `alpha` for a one-sided
alternative and `alpha / 2` for a two-sided alternative. Only the
starting value of the search is affected, so the returned sample size is
exact in either case.

The search stops where the exact power crosses the target: the returned
size of group 2 attains the target power and the size one unit smaller
does not, unless the returned size is 1. The exact power is not monotone
in the sample size, so a smaller size further from the normal
approximation may also attain the target power, and a larger size may
fall short of it.

The two other searches do not depend on where the normal approximation
starts. `search = 'smallest'` scans the size of group 2 upwards from 1
and returns the first size that attains the target power. It visits
every smaller size, which takes longer for the unconditional tests.
`search = 'stable'` scans downwards from the limit given by
`search.limit` and returns the smallest size from which every size up to
the limit attains the target power, so that adding patients within this
range does not take the power below the target. The search stops with an
error if the limit itself falls short of the target power. The three
searches can return different sizes because the exact power is not
monotone in the sample size.

Under `method = 'standard'`, `'null.variance'` or
`'alternative.variance'` the steps above are replaced by the closed-form
normal approximation and the rounding rule selected by `rounding`. The
`Power` column still reports the exact power of `Test` at the resulting
sample size.

With a non-zero `margin`, the difference `p1 - p2` in the normal
approximation is replaced by its distance `p1 - p2 + margin` from the
boundary of the null hypothesis, and the variance under the null
hypothesis is evaluated at the large sample values of the restricted
maximum likelihood estimates of Farrington and Manning (1990). The
method `'standard'` then gives formula (4) of Farrington and Manning
(1990), which is formula (2) of Friede et al. (2007), and
`'alternative.variance'` gives the formula of Blackwelder (1982), which
is formula (1) of Friede et al. (2007). The exact search starts from the
formula of Farrington and Manning.

With `margin.scale = 'RR'` and the ratio `R0 = margin`, the difference
in the normal approximation is the distance `p1 - R0 p2`, or
`R0 p2 - p1` for `alternative = 'less'`, and the variance of the size of
group 2 is `p1 (1 - p1) / r + R0^2 p2 (1 - p2)`, evaluated at the large
sample values of the restricted maximum likelihood estimates under
`p1 = R0 p2` for the significance term. The method `'standard'` then
gives formula (8) of Farrington and Manning (1990),
`'alternative.variance'` gives the formula of their Method 1, and the
exact search starts from formula (8).

## References

Friede T, Kieser M (2004). Sample size recalculation for binary data in
internal pilot study designs. *Pharmaceutical Statistics*, 3(4),
269-279.

Kieser M (2020). *Methods and Applications of Sample Size Calculation
and Recalculation in Clinical Trials*. Springer, Cham.

Blackwelder WC (1982). "Proving the null hypothesis" in clinical trials.
*Controlled Clinical Trials*, 3(4), 345-353.

Farrington CP, Manning G (1990). Test statistics and sample size
formulae for comparative binomial trials with null hypothesis of
non-zero risk difference or non-unity relative risk. *Statistics in
Medicine*, 9(12), 1447-1454.

Friede T, Mitchell C, Mueller-Velten G (2007). Blinded sample size
reestimation in non-inferiority trials with binary endpoints.
*Biometrical Journal*, 49(6), 903-916.

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
# One-sided chi-squared test
BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
                 tar.power = 0.8, Test = 'Chisq')
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Chisq
#>   Alternative      : greater
#>   Response rates   : p1 = 0.4, p2 = 0.2
#>   Allocation ratio : 1 to 1
#>   Alpha            : 0.025
#>   Target power     : 0.8
#>   Exact search     : crossing
#> 
#>   Required sample size: N1 = 80, N2 = 80, total N = 160
#>   Attained power      : 0.8009

# \donttest{
# Two-sided Fisher exact test
BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 2, alpha = 0.05,
                 tar.power = 0.9, Test = 'Fisher', alternative = 'two.sided')
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Fisher
#>   Alternative      : two.sided
#>   Response rates   : p1 = 0.5, p2 = 0.2
#>   Allocation ratio : 2 to 1
#>   Alpha            : 0.05
#>   Target power     : 0.9
#>   Exact search     : crossing
#> 
#>   Required sample size: N1 = 78, N2 = 39, total N = 117
#>   Attained power      : 0.9013

# Normal approximation with the variance under the null hypothesis, two-sided test
# at level 0.05, as in Friede and Kieser (2004)
BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.05, tar.power = 0.8,
                 Test = 'Chisq', alternative = 'two.sided',
                 method = 'null.variance')
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Chisq
#>   Alternative      : two.sided
#>   Response rates   : p1 = 0.4, p2 = 0.2
#>   Allocation ratio : 1 to 1
#>   Alpha            : 0.05
#>   Target power     : 0.8
#> 
#>   Required sample size: N1 = 83, N2 = 83, total N = 166
#>   Attained power      : 0.8137

# Non-inferiority with a margin of 0.1, formula (2) of Friede et al. (2007)
BinarySampleSize(p1 = 0.7, p2 = 0.7, r = 1, alpha = 0.025, tar.power = 0.8,
                 Test = 'Farrington-Manning', method = 'standard',
                 rounding = 'friede-kieser', margin = 0.1)
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Farrington-Manning
#>   Alternative      : greater
#>   Response rates   : p1 = 0.7, p2 = 0.7
#>   Allocation ratio : 1 to 1
#>   Margin           : 0.1
#>   Alpha            : 0.025
#>   Target power     : 0.8
#> 
#>   Required sample size: N1 = 329, N2 = 329, total N = 658
#>   Attained power      : 0.8006

# Non-inferiority on the scale of the risk ratio: the null hypothesis p1 / p2 <= 0.8
# with equal response probabilities of 0.6
BinarySampleSize(p1 = 0.6, p2 = 0.6, r = 1, alpha = 0.025, tar.power = 0.8,
                 Test = 'Farrington-Manning', margin = 0.8, margin.scale = 'RR')
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Farrington-Manning
#>   Alternative      : greater
#>   Response rates   : p1 = 0.6, p2 = 0.6
#>   Allocation ratio : 1 to 1
#>   Margin           : 0.8 (risk ratio)
#>   Alpha            : 0.025
#>   Target power     : 0.8
#>   Exact search     : crossing
#> 
#>   Required sample size: N1 = 215, N2 = 215, total N = 430
#>   Attained power      : 0.802
# }
```
