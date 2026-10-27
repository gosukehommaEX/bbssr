# Blinded Sample Size Re-Estimation from Paired Interim Data

Re-estimates the number of pairs of a matched-pairs trial from the
number of discordant pairs observed at an interim analysis, without the
direction in which they are discordant.

## Usage

``` r
PairedBSSR(
  n1,
  D1,
  Delta.A,
  alpha,
  tar.power,
  Test,
  restricted = FALSE,
  N = NULL,
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
  margin = 0,
  ref.pvalue = FALSE
)
```

## Arguments

- n1:

  Number of pairs at the interim analysis

- D1:

  Number of discordant pairs among them

- Delta.A:

  Assumed difference `p10 - p01` used in the re-estimation

- alpha:

  Level of significance for the alternative specified by `alternative`

- tar.power:

  Target power of the re-estimation

- Test:

  Type of statistical test. Options: `'McNemar'`, `'McNemar-cc'`,
  `'McNemar-exact'`, `'McNemar-midP'`, `'McNemar-uncond'` or
  `'Paired-Score'`

- restricted:

  Logical. If `TRUE`, the final number of pairs is at least `N`. Default
  is `FALSE`

- N:

  Planned number of pairs. Required when `restricted = TRUE`, otherwise
  optional

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

  Level used in the re-estimation. Default is `alpha`

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

An object of class `bbssr_pbssr`, a data frame with one row containing
the interim data `n1` and `D1`, the blinded estimate `psi.hat` of the
probability of a discordant pair, the value `psi.plan` used in the
re-estimation, the planning probabilities `p10.plan` and `p01.plan`, the
re-estimated number of pairs `N.re` before the bounds are applied (at
most `N.max` with `ss.method = 'exact'`, whose search does not go beyond
it, and within the bounds with `ss.method = 'conditional.power'`), the
number `n.stage2` of pairs still to be enrolled, the final number of
pairs `N.final`, the power `Power` of the final analysis with `N.final`
pairs under the planning probabilities, and its conditional power `CP`
given `D1`, computed with `ss.Test` at the level `ss.alpha`. The
attributes hold the design settings.

## Details

The re-estimation is that of
[`PairedPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedPowerBSSR.md):
the blinded estimate `D1 / n1` of the probability of a discordant pair,
raised to `|Delta.A|` when it is smaller, and the assumed difference
`Delta.A` give the planning probabilities, from which the number of
pairs follows by the rule `ss.method`. The conditional power is
`sum over D2 of b(D2; N.final - n1, psi.plan) g(D1 + D2)`, the
probability of rejection given `D1` under the planning probabilities,
formula (4) of Schepers et al. (2026) for the exact conditional McNemar
test.

## References

Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026).
Blinded sample size review for McNemar's test based on primary and
surrogate endpoints. arXiv:2608.17784.

## Examples

``` r
# 25 interim pairs of which 5 are discordant, assumed difference 0.15
PairedBSSR(n1 = 25, D1 = 5, Delta.A = 0.15, alpha = 0.05, tar.power = 0.8,
           Test = 'McNemar-exact', alternative = 'two.sided',
           ss.method = 'conditional.power', N.max = 150)
#> Blinded sample size re-estimation from paired interim data
#> 
#>   Test             : McNemar-exact
#>   Alternative      : two.sided
#>   Re-estimation    : conditional power (unrestricted)
#>   Assumed effect   : 0.15
#>   Alpha            : 0.05, target power 0.8
#> 
#> Interim data
#>   Pairs            : n1 = 25, discordant D1 = 5
#>   Discordance      : estimate 0.2, used 0.2
#>   Planning values  : p10 = 0.175, p01 = 0.025
#> 
#> Re-estimation
#>   Re-estimated N   : 70
#>   Still to enrol   : 45 pairs
#>   Final N          : 70
#>   Power at final N : 0.7911
#>   Conditional power: 0.8019
```
