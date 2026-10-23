# Power of Blinded Sample Size Re-Estimation for Paired Binary Data

Calculates the exact power and the distribution of the final number of
pairs of a matched-pairs design with blinded sample size re-estimation,
together with the power of the fixed design with the planned number of
pairs. The number of pairs is re-estimated at an interim analysis from
the number of discordant pairs, without the direction in which they are
discordant, and the final analysis uses all pairs.

## Usage

``` r
PairedPowerBSSR(
  psi,
  Delta.T,
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
  margin = 0,
  ref.pvalue = FALSE
)
```

## Arguments

- psi:

  True probability `p10 + p01` of a discordant pair. A vector is allowed

- Delta.T:

  True difference `p10 - p01`, a single value

- Delta.A:

  Assumed difference `p10 - p01` used in the re-estimation

- N:

  Planned number of pairs, which is the number of pairs of the fixed
  design

- omega:

  Fraction of `N` at the interim analysis, rounded up. Supply either
  `omega` or `n.interim`

- alpha:

  Level of significance for the alternative specified by `alternative`

- tar.power:

  Target power of the re-estimation

- Test:

  Type of statistical test. Options: `'McNemar'`, `'McNemar-exact'`,
  `'McNemar-midP'`, `'McNemar-uncond'` or `'Paired-Score'`

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

  Level used in the re-estimation. Default is `alpha`

- search:

  Search of the exact re-estimation, see
  [`PairedSampleSize`](https://gosukehommaex.github.io/bbssr/reference/PairedSampleSize.md)

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

An object of class `bbssr_ppowerbssr`, a data frame with one row for
each element of `psi` that gives probabilities `p10` and `p01` in
`[0, 1]`, containing:

- psi, p10, p01:

  True probabilities of a discordant pair and of its two kinds

- power.BSSR:

  Power of the design with re-estimation

- power.TRAD:

  Power of the fixed design with `N` pairs

- E.N:

  Expected final number of pairs

The attribute `reestimation` holds the final number of pairs for every
interim number of discordant pairs, `N.dist` the distribution of the
final number of pairs in each scenario, and the other attributes the
design settings.

## Details

Among the `n1` pairs of the interim analysis, the number `D1` of
discordant pairs is binomial with `n1` trials and probability `psi`, and
the blinded estimate of `psi` is `D1 / n1`. The assumed difference
`Delta.A` is kept, and the estimate is raised to `|Delta.A|` when it is
smaller, since `psi >= |p10 - p01|`. The planning probabilities are then
`p10 = (psi.plan + Delta.A) / 2` and `p01 = (psi.plan - Delta.A) / 2`,
as in Xu and Yu (2013) and Schepers et al. (2026).

The final number of pairs depends on the interim data only through `D1`.
Given the numbers of discordant pairs of the two stages, the number of
pairs in which only treatment 1 responds is binomial with the total
number of discordant pairs and probability `p10 / psi`, whatever the
re-estimation did. The power is therefore
`sum over D1 and D2 of b(D1; n1, psi) b(D2; N(D1) - n1, psi) g(D1 + D2)`,
where `g(D)` is the probability that the final analysis with `N(D1)`
pairs rejects given `D` discordant pairs. It is computed exactly by
enumeration.

With `ss.method = 'conditional.power'` the final number of pairs is the
smallest number `N` between the bounds at which
`sum over D2 of b(D2; N - n1, psi.plan) g_N(D1 + D2)`, computed with the
planning probabilities and the test `ss.Test` at the level `ss.alpha`,
attains `tar.power`. With the exact conditional McNemar test this is the
conditional power (4) of Schepers et al. (2026), which they evaluate
with `p* = 0.5 + Delta / (2 psi)` at the blinded estimate of `psi`. When
no number up to `N.max` attains the target, `N.max` is used, a choice
the article leaves open.

The exact conditional McNemar test keeps its level whatever blinded rule
is used, since the number of pairs in which only treatment 1 responds
stays binomial given the total number of discordant pairs (Schepers et
al., 2026). The other tests can exceed their nominal level after
re-estimation, see
[`PairedTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/PairedTypeIErrorBSSR.md).

## References

Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026).
Blinded sample size review for McNemar's test based on primary and
surrogate endpoints. arXiv:2608.17784.

Xu J, Yu M (2013). Sample size determination and re-estimation for
matched pair designs with multiple binary endpoints. *Biometrical
Journal*, 55(3), 430-443.

## Examples

``` r
res <- PairedPowerBSSR(psi = c(0.15, 0.2, 0.3), Delta.T = 0.15, Delta.A = 0.15, N = 72,
                       n.interim = 30, alpha = 0.05, tar.power = 0.8,
                       Test = 'McNemar-exact', alternative = 'two.sided',
                       ss.method = 'standard')
print(res)
#> Blinded sample size re-estimation for paired binary data
#> 
#>   Test             : McNemar-exact
#>   Alternative      : two.sided
#>   Re-estimation    : normal approximation (unrestricted)
#>   Pairs            : N = 72 planned, 30 at the interim analysis
#>   Effect           : assumed 0.15, true 0.15
#>   Alpha            : 0.05, target power 0.8
#> 
#>   psi   p10   p01 power.BSSR power.TRAD    E.N
#>  0.15 0.150 0.000     0.7851     0.9684  59.33
#>  0.20 0.175 0.025     0.7174     0.8056  71.22
#>  0.30 0.225 0.075     0.7190     0.5824 102.98

# \donttest{
# Conditional power rule of Schepers et al. (2026) with at most 150 pairs
res <- PairedPowerBSSR(psi = c(0.15, 0.2, 0.3), Delta.T = 0.15, Delta.A = 0.15, N = 72,
                       n.interim = 30, alpha = 0.05, tar.power = 0.8,
                       Test = 'McNemar-exact', alternative = 'two.sided',
                       ss.method = 'conditional.power', N.max = 150)
print(res)
#> Blinded sample size re-estimation for paired binary data
#> 
#>   Test             : McNemar-exact
#>   Alternative      : two.sided
#>   Re-estimation    : conditional power (unrestricted)
#>   Pairs            : N = 72 planned, 30 at the interim analysis
#>   Effect           : assumed 0.15, true 0.15
#>   Alpha            : 0.05, target power 0.8
#> 
#>   psi   p10   p01 power.BSSR power.TRAD    E.N
#>  0.15 0.150 0.000     0.8946     0.9684  63.98
#>  0.20 0.175 0.025     0.7415     0.8056  74.92
#>  0.30 0.225 0.075     0.7428     0.5824 108.17
# }
```
