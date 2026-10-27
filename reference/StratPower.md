# Exact Power of a Stratified Test

Calculates the exact power of a test of
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md)
for a two-arm trial with a binary endpoint, with fixed group sizes in
each stratum or, with two strata, with stratum and group sizes that are
random.

## Usage

``` r
StratPower(
  p1,
  p2,
  N1 = NULL,
  N2 = NULL,
  alpha,
  Test,
  alternative = c("greater", "less", "two.sided"),
  tsmethod = c("minlike", "central", "blaker"),
  alloc = c("fixed", "stratum", "group", "both"),
  N = NULL,
  frac = NULL,
  r = 1,
  block = 4
)
```

## Arguments

- p1:

  Response probabilities of group 1 in the strata: a vector with one
  element per stratum, or a matrix with one row per scenario and one
  column per stratum

- p2:

  Response probabilities of group 2 in the strata, of the same form as
  `p1`

- N1:

  Sizes of group 1 in the strata, one element per stratum, with
  `alloc = 'fixed'`

- N2:

  Sizes of group 2 in the strata, one element per stratum, with
  `alloc = 'fixed'`

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

- alloc:

  Allocation of the patients to the strata and the groups. `'fixed'`
  (default) uses the group sizes `N1` and `N2`. The others need two
  strata and make the sizes random, see Details: `'stratum'` (stratum
  sizes binomial, permuted blocks within the strata), `'group'` (stratum
  sizes fixed, group sizes binomial) and `'both'` (stratum and group
  sizes binomial)

- N:

  Total sample size, with a random allocation

- frac:

  Fractions of the two strata, with a random allocation. They are
  rescaled to sum to 1

- r:

  Allocation ratio of group 1 to group 2 within the strata, a single
  value or one value per stratum, with a random allocation. Default is 1

- block:

  Size of the permuted blocks of `alloc = 'stratum'`. Each block holds
  `block r / (1 + r)` patients of group 1, which must be a whole number.
  Default is 4

## Value

An object of class `bbssr_spower`, a data frame with one row per
scenario containing the response probabilities `p1.S1, p1.S2, ...` and
`p2.S1, p2.S2, ...` and the exact power `Power`. With a random
allocation the attribute `skipped` is the total probability of the
configurations left out of the sum, see Details, and `n.strata` the
stratum sizes of `alloc = 'group'`.

## Details

The numbers of responders of the two groups in each stratum are
independent binomial variables. The power is the probability of the
rejection region of
[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md),
computed exactly by enumerating the numbers of responders of the strata
and carrying the distribution of the partial sum of `S` over the strata.
The response probabilities may differ between the strata in any way, so
the power is available also when the odds ratio or the risk difference
is not common to the strata.

With `alloc` other than `'fixed'` the design has two strata and a total
of `N` patients, and the sizes of the strata and of the groups within
them are random. With `'group'` the stratum sizes are fixed at
`n1 = floor(N frac[1])` and `n2 = N - n1`, and each patient of stratum
`j` is in group 1 with the probability `r[j] / (1 + r[j])`, so that the
size of group 1 is binomial (Jung, 2014, Section 3.2, where the size of
stratum 1 is the round-off of `N frac[1]`; the integer part used here
reproduces all but one of the sample sizes of the column 'n fixed' of
Table 4 of the article). With `'both'` the size of stratum 1 is binomial
as well, with the probability `frac[1]` (Jung, 2014, Section 3.1). With
`'stratum'` the stratum sizes are binomial and the groups are allocated
by permuted blocks of size `block` within each stratum: every complete
block holds `block r[j] / (1 + r[j])` patients of group 1, and the
number of patients of group 1 among the patients of the last, incomplete
block is hypergeometric. The power is the sum over the stratum and group
sizes of their probability times the power of the fixed design they
give. Configurations with a probability below 1e-14 are left out, and
their total probability, which bounds the error of the power, is
returned in the attribute `skipped`. A group without patients in a
stratum leaves that stratum without variation of `S`, and a stratum
without patients leaves a design with one stratum. The number of
configurations grows with `N` most quickly for `'both'`: at `N = 100`
the power of one scenario takes a fraction of a second with `'stratum'`,
about a second with `'group'` and tens of seconds with `'both'`.

## References

Jung SH (2014). Stratified Fisher's exact test and its sample size
calculation. *Biometrical Journal*, 56(1), 129-140.

## See also

[`StratRR`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md),
[`StratSampleSize`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md)

## Author

Gosuke Homma (<my.name.is.gosuke@gmail.com>)

## Examples

``` r
StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), N1 = c(40, 40), N2 = c(40, 40),
           alpha = 0.025, Test = 'CMH')
#> Exact power of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Group sizes      : N1 = (40, 40), N2 = (40, 40)
#>   Alpha            : 0.025
#> 
#>  p1.S1 p1.S2 p2.S1 p2.S2  Power
#>    0.5   0.7   0.3   0.5 0.7322

# Several scenarios, and the stratified exact test
StratPower(p1 = rbind(c(0.5, 0.7), c(0.45, 0.65)), p2 = rbind(c(0.3, 0.5), c(0.3, 0.5)),
           N1 = c(40, 40), N2 = c(40, 40), alpha = 0.025, Test = 'Strat-Fisher')
#> Exact power of a stratified test
#> 
#>   Test             : Strat-Fisher
#>   Alternative      : greater
#>   Group sizes      : N1 = (40, 40), N2 = (40, 40)
#>   Alpha            : 0.025
#> 
#>  p1.S1 p1.S2 p2.S1 p2.S2  Power
#>   0.50  0.70   0.3   0.5 0.6792
#>   0.45  0.65   0.3   0.5 0.4249

# Random stratum sizes, with permuted blocks of 4 within the strata
StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), alpha = 0.025, Test = 'CMH',
           alloc = 'stratum', N = 160, frac = c(0.5, 0.5))
#> Exact power of a stratified test
#> 
#>   Test             : CMH
#>   Alternative      : greater
#>   Allocation       : random stratum sizes, permuted blocks of 4 within the strata
#>   Total size       : N = 160, fractions (0.5, 0.5), ratio 1
#>   Left out         : configurations of total probability 8.4e-14
#>   Alpha            : 0.025
#> 
#>  p1.S1 p1.S2 p2.S1 p2.S2  Power
#>    0.5   0.7   0.3   0.5 0.7332
```
