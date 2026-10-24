# Rejection Region of a Stratified Test

Returns the critical values of a test for a two-arm trial with a binary
endpoint stratified into two or more strata, for every combination of
the numbers of responders of the strata. The tests are based on `S`, the
number of responders of group 1 summed over the strata. Given the
numbers of responders `z = (z_1, ..., z_J)` of the strata, each test
rejects the null hypothesis when `S` is at most a lower or at least an
upper critical value. The unstratified tests, which ignore the strata,
are included for comparison.

## Usage

``` r
StratRR(
  N1,
  N2,
  alpha,
  Test,
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker")
)
```

## Arguments

- N1:

  Sizes of group 1 in the strata, one element per stratum

- N2:

  Sizes of group 2 in the strata, one element per stratum

- alpha:

  Level of significance, for the alternative specified

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

## Value

An object of class `bbssr_srr`, a data frame with one row for each
combination of the numbers of responders of the strata, the first
stratum varying fastest, containing:

- z.S1, z.S2, ...:

  Numbers of responders of the strata

- S.min, S.max:

  Smallest and largest possible value of `S`

- lower:

  Lower critical value: the null hypothesis is rejected when
  `S <= lower`. `NA` when no value of `S` is rejected in the lower tail

- upper:

  Upper critical value: the null hypothesis is rejected when
  `S >= upper`. `NA` when no value of `S` is rejected in the upper tail

## Details

In stratum `j`, with `N_j = N1[j] + N2[j]` patients and `z_j`
responders, the null expectation of the number of responders of group 1
is `N1[j] z_j / N_j`, and `E` is the sum of these over the strata. The
Cochran-Mantel-Haenszel statistic is `(S - E) / sqrt(V)` with the
hypergeometric variance
`V = sum_j N1[j] N2[j] z_j (N_j - z_j) / (N_j^2 (N_j - 1))`, and
Cochran's statistic uses the binomial variance
`sum_j N1[j] N2[j] z_j (N_j - z_j) / N_j^3` instead. The continuity
correction reduces `|S - E|` by 0.5 when it is at least 0.5, as in
[`mantelhaen.test`](https://rdrr.io/r/stats/mantelhaen.test.html). Nam
(1992) reduces `S - E` by 0.5 in every outcome, which makes a difference
only where `|S - E| < 0.5`. There `S` is the integer next to `E`, and
since the hypergeometric variance is at least `f (1 - f)`, with `f` the
fractional part of `E`, and the binomial variance at least half of it,
the one-sided p-value of `'Cochran-cc'` exceeds 0.078. The two
corrections therefore reject the same outcomes at one-sided levels up to
0.078. The statistics are referred to the standard normal distribution,
with the absolute value for a two-sided alternative, and an outcome
whose variance is 0 is not rejected. The stratified exact conditional
test uses the distribution of `S` given `z`, the convolution of the
hypergeometric distributions of the strata. The two-sided p-value
`'central'` is twice the smaller tail probability (Jung, 2014), and
`'minlike'` and `'blaker'` order the values of `S` by their probability
and by their smaller tail probability, with the ties of
[`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md).
As in the other functions of the package, a p-value is rejected when it
is below `alpha` by more than about 1.49e-8.

The unstratified tests apply the corresponding test of a single 2 x 2
table to the table of totals, with `sum_j N1[j]` and `sum_j N2[j]`
patients and `sum_j z_j` responders, so their critical values depend on
`z` only through this total. Pearson's statistic of that table is
Cochran's statistic of one stratum, as in `'Chisq'` of
[`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md),
and the continuity correction is that of `'Cochran-cc'`; with the
correction applied to every outcome, as in the \\z_u\\ test of Nam
(1992, Section 2.2), the same outcomes are rejected at one-sided levels
up to 0.078. An unstratified test is valid when the distribution of the
strata is the same in the two groups, and its level can be distorted
otherwise (Nam, 1992; Jung, 2014); its power and type I error rate are
computed here with the response probabilities of each stratum, so these
consequences can be evaluated.

The function stores the critical values for reuse during the session,
like the rejection regions of the other designs.

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
[`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
rr <- StratRR(N1 = c(10, 12), N2 = c(10, 12), alpha = 0.025, Test = 'CMH')
print(rr)
#> Rejection region of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Group sizes      : N1 = (10, 12), N2 = (10, 12)
#>   Alpha            : 0.025
#>   Combinations     : 525, of which 485 have a rejected value of S
#> 
#>  z.S1 z.S2 S.min S.max lower upper
#>     4    0     0     4    NA     4
#>     5    0     0     5    NA     5
#>     6    0     0     6    NA     6
#>     7    0     0     7    NA     6
#>     8    0     0     8    NA     7
#>     9    0     0     9    NA     7
#>    10    0     0    10    NA     8
#>    11    0     1    10    NA     8
#>    12    0     2    10    NA     9
#>    13    0     3    10    NA     9
#>    14    0     4    10    NA    10
#>    15    0     5    10    NA    10
#>    16    0     6    10    NA    10
#>     3    1     0     4    NA     4
#>     4    1     0     5    NA     5
#>     5    1     0     6    NA     6
#>     6    1     0     7    NA     6
#>     7    1     0     8    NA     7
#>     8    1     0     9    NA     7
#>     9    1     0    10    NA     8
#>   ... 465 more rows
```
