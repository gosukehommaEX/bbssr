# Sample Size for Paired Binary Data

Calculates the number of pairs required for a test comparing two
treatments on a binary endpoint in a matched-pairs design to attain a
target power. The number can be obtained from the exact power of the
selected test or from the normal approximation.

## Usage

``` r
PairedSampleSize(
  p10,
  p01,
  alpha,
  tar.power,
  Test,
  alternative = c("greater", "less", "two.sided", "equivalence"),
  n.grid = 100,
  bb.gamma = 0,
  method = c("exact", "standard"),
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50),
  margin = 0,
  ref.pvalue = FALSE
)
```

## Arguments

- p10:

  Probability that only treatment 1 responds in a pair

- p01:

  Probability that only treatment 2 responds in a pair

- alpha:

  Level of significance for the alternative specified by `alternative`

- tar.power:

  Target power

- Test:

  Type of statistical test. Options: `'McNemar'`, `'McNemar-exact'`,
  `'McNemar-midP'`, `'McNemar-uncond'` or `'Paired-Score'`

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

- method:

  How the number of pairs is obtained. `'exact'` (default) steps from
  the normal approximation to the number of pairs at which the exact
  power of `Test` attains the target, see Details. `'standard'` rounds
  up the normal approximation, formula (9) of Nam (1997), which is
  formula (3) of Connor (1987) when the margin is 0

- search:

  How the search of `method = 'exact'` chooses the number of pairs.
  `'crossing'` (default) steps from the normal approximation to a number
  that attains the target power while the number one smaller does not.
  `'smallest'` returns the smallest number that attains the target
  power. `'stable'` returns the smallest number from which every number
  up to a limit attains it. Ignored by `method = 'standard'`

- search.limit:

  Two numbers `a` and `b` giving the limit
  `max(ceiling(a n0), ceiling(n0 + b))` of the number of pairs examined
  by `search = 'stable'`, where `n0` is the normal approximation rounded
  up. `a` must be at least 1 and `b` non-negative. Default is `c(2, 50)`

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

An object of class `bbssr_psamplesize`, a data frame with one row
containing:

- p10, p01:

  Probabilities of the two kinds of discordant pairs

- alpha:

  Level of significance

- tar.power:

  Target power

- Test:

  Name of the statistical test

- alternative:

  Direction of the alternative hypothesis

- Power:

  Exact power of `Test` with the required number of pairs

- N:

  Required number of pairs

The attributes hold the remaining settings, among them the unrounded
normal approximation `n.raw` and the limit `search.limit` of the search
`'stable'` (otherwise `NA`).

## Details

The normal approximation is formula (9) of Nam (1997),
`n = (z.a sqrt(v0) + z.b sqrt(v1))^2 / (d1 - d0)^2`, where
`d1 = p10 - p01` is the true difference, `d0 = -margin` the difference
under the null hypothesis, `v1 = p10 + p01 - d1^2` the variance under
the alternative, and `v0 = 2 p01.bar + d0 - d0^2` the variance under the
null hypothesis evaluated at the large sample value `p01.bar` of the
restricted maximum likelihood estimate of `p01`, formula (5) of that
article. `z.a` is the upper `alpha` quantile of the standard normal
distribution, the upper `alpha / 2` quantile for a two-sided
alternative, and `z.b` the `tar.power` quantile. With `margin = 0` the
formula is formula (3) of Connor (1987). For `alternative = 'less'` the
probabilities `p10` and `p01` are exchanged. For
`alternative = 'equivalence'` the approximate power of each one-sided
test is obtained in the same way for its own margin, at the level
`alpha`, and the normal approximation is the smallest number of pairs at
which the sum of the two approximate powers minus 1 attains the target
power. The sum minus 1 is a lower bound of the approximate probability
that both tests reject.

The normal approximation is the starting point of the exact search for
every test. The exact power is not monotone in the number of pairs, so
the three searches can return different numbers, see
[`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).

## References

Connor RJ (1987). Sample size for testing differences in proportions for
the paired-sample design. *Biometrics*, 43(1), 207-211.

Nam J (1997). Establishing equivalence of two treatments and sample size
requirements in matched-pairs design. *Biometrics*, 53(4), 1422-1430.

## Examples

``` r
# Example of Connor (1987, p. 210): theta1 = 0.8, theta2 = 0.7 and psi = 0.3
PairedSampleSize(p10 = 0.2, p01 = 0.1, alpha = 0.05, tar.power = 0.9,
                 Test = 'McNemar', method = 'standard')
#> Number of pairs for paired binary data
#> 
#>   Test           : McNemar
#>   Alternative    : greater
#>   Discordance    : p10 = 0.2, p01 = 0.1
#>   Alpha          : 0.05
#>   Target power   : 0.9
#>   Method         : normal approximation
#> 
#>   Required number of pairs: N = 254
#>   Attained power          : 0.9062

# \donttest{
# Exact number of pairs of the exact conditional test
PairedSampleSize(p10 = 0.2, p01 = 0.1, alpha = 0.05, tar.power = 0.9,
                 Test = 'McNemar-exact')
#> Number of pairs for paired binary data
#> 
#>   Test           : McNemar-exact
#>   Alternative    : greater
#>   Discordance    : p10 = 0.2, p01 = 0.1
#>   Alpha          : 0.05
#>   Target power   : 0.9
#>   Exact search   : crossing
#> 
#>   Required number of pairs: N = 269
#>   Attained power          : 0.9
# }
```
