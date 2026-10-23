# Sample Size of a Stratified Test

Calculates the total sample size of a two-arm trial with a binary
endpoint stratified into two or more strata, for a test of
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md),
either exactly or by the normal approximation of Woolson et al. (1986)
and Nam (1992).

## Usage

``` r
StratSampleSize(
  p1,
  p2,
  frac,
  r = 1,
  alpha,
  tar.power,
  Test,
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker"),
  method = c("exact", "standard"),
  search = c("crossing", "smallest", "stable"),
  search.limit = c(2, 50)
)
```

## Arguments

- p1:

  Response probabilities of group 1 in the strata, one element per
  stratum

- p2:

  Response probabilities of group 2 in the strata, one element per
  stratum

- frac:

  Fractions of the strata, one positive element per stratum. They are
  rescaled to sum to 1, so they may also be given as numbers of patients

- r:

  Allocation ratio of group 1 to group 2 within the strata, a single
  value or one value per stratum. Default is 1

- alpha:

  Level of significance, for the alternative specified

- tar.power:

  Target power

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

- method:

  Calculation of the sample size. `'exact'` (default) searches for the
  total sample size at which the exact power of
  [`StratPower`](https://gosukehommaex.github.io/bbssr/reference/StratPower.md)
  attains `tar.power`. `'standard'` uses the normal approximation, see
  Details

- search:

  How the exact calculation chooses among the sizes at which the power
  fluctuates around `tar.power`, as in
  [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).
  Options: `'crossing'` (default), `'smallest'` or `'stable'`

- search.limit:

  Limit of the search `'stable'`: the nominal sizes are examined
  downwards from
  `max(ceiling(search.limit[1] n), ceiling(n + search.limit[2]))`, where
  `n` is the nominal size of the normal approximation. Default is
  `c(2, 50)`

## Value

An object of class `bbssr_ssamplesize`, a list containing the arguments
`p1`, `p2`, `frac`, `r`, `alpha`, `tar.power`, `Test` and `alternative`,
the group sizes `N1` and `N2` of the strata, their total `N`, and the
exact power `Power` of the design. The attribute `n.raw` holds the
unrounded total of the normal approximation, `n` the nominal total of
the design and `search.limit` the limit of the nominal total examined by
`search = 'stable'` (otherwise `NA`).

## Details

A design is described by the fractions `frac` of the strata and the
allocation ratios `r`. For a nominal total sample size `n`, the size of
group 2 in stratum `j` is `ceiling(frac[j] n / (1 + r[j]))` and that of
group 1 is `ceiling(r[j] N2[j])`, so the total `N` is at least `n`. The
exact calculation searches over `n` within this family of designs.

The normal approximation is that of Cochran's test. With
`s[j] = r[j] / (1 + r[j])`, `w[j] = frac[j] s[j] (1 - s[j])`,
`D = sum_j w[j] (p1[j] - p2[j])`,
`v0 = sum_j w[j] pbar[j] (1 - pbar[j])` with
`pbar[j] = s[j] p1[j] + (1 - s[j]) p2[j]`, and
`v1 = sum_j w[j] ((1 - s[j]) p1[j] (1 - p1[j]) + s[j] p2[j] (1 - p2[j]))`,
the total sample size is `N* = (z.a sqrt(v0) + z.b sqrt(v1))^2 / D^2`,
formula (2.8) of Woolson et al. (1986), where `z.a` is the upper `alpha`
quantile (`alpha / 2` for a two-sided alternative) and `z.b` the
`tar.power` quantile of the standard normal distribution. For the tests
with a continuity correction it is
`N = (1 + sqrt(1 + 2 / (D N*)))^2 N* / 4`, formula (2.3) of Nam (1992).
The same approximation is used for every test, since the
Cochran-Mantel-Haenszel statistic differs from Cochran's only by the
factors `(N_j - 1) / N_j` of the variance, and it serves as the starting
value of the exact search. The nominal total is `ceiling(N*)`, or
`ceiling(N)`; the published tables of Woolson et al. (1986) and Nam
(1992) compute `N*` and `N` with the quantiles rounded to three
decimals, such as 1.645 for the level 0.05 and 0.841 or 0.842 for the
power 0.8, and round them to the nearest integer instead.

The exact power is not monotone in the sample size, so the result
depends on `search`, as described in
[`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).

## References

Woolson RF, Bean JA, Rojas PB (1986). Sample size for case-control
studies using Cochran's statistic. *Biometrics*, 42(4), 927-932.

Nam J (1992). Sample size determination for case-control studies and the
comparison of stratified and unstratified analyses. *Biometrics*, 48(2),
389-395.

Jung SH (2014). Stratified Fisher's exact test and its sample size
calculation. *Biometrical Journal*, 56(1), 129-140.

## See also

[`StratPower`](https://gosukehommaex.github.io/bbssr/reference/StratPower.md),
[`StratPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
# Normal approximation and exact sample size with two strata
StratSampleSize(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.5, 0.5),
                alpha = 0.025, tar.power = 0.8, Test = 'CMH', method = 'standard')
#> Sample size of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Method           : normal approximation
#>   Alpha            : 0.025, target power 0.8
#>   Normal approx.   : 186 (unrounded total)
#>   Total size       : N = 188, power 0.7985
#> 
#>  Stratum frac r  p1  p2 N1 N2
#>        1  0.5 1 0.5 0.3 47 47
#>        2  0.5 1 0.7 0.5 47 47
StratSampleSize(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.5, 0.5),
                alpha = 0.025, tar.power = 0.8, Test = 'CMH')
#> Sample size of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Method           : exact (search crossing)
#>   Alpha            : 0.025, target power 0.8
#>   Normal approx.   : 186 (unrounded total)
#>   Total size       : N = 192, power 0.8085
#> 
#>  Stratum frac r  p1  p2 N1 N2
#>        1  0.5 1 0.5 0.3 48 48
#>        2  0.5 1 0.7 0.5 48 48
```
