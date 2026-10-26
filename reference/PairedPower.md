# Power for Paired Binary Data

Calculates the exact power of a test comparing two treatments on a
binary endpoint in a matched-pairs design, for given probabilities of
the two kinds of discordant pairs.

## Usage

``` r
PairedPower(
  p10,
  p01,
  N,
  alpha,
  Test,
  alternative = c("greater", "less", "two.sided", "equivalence"),
  n.grid = 100,
  bb.gamma = 0,
  margin = 0,
  ref.pvalue = FALSE
)
```

## Arguments

- p10:

  Probability that only treatment 1 responds in a pair. A vector is
  allowed

- p01:

  Probability that only treatment 2 responds in a pair, of the same
  length as `p10`. A single value of either is used with every value of
  the other

- N:

  Number of pairs

- alpha:

  Level of significance for the alternative specified by `alternative`

- Test:

  Type of statistical test. Options: `'McNemar'`, `'McNemar-cc'`,
  `'McNemar-exact'`, `'McNemar-midP'`, `'McNemar-uncond'` or
  `'Paired-Score'`

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

An object of class `bbssr_ppower`, a data frame with one row for each
element of `p10`, containing:

- p10, p01:

  Probabilities of the two kinds of discordant pairs

- N:

  Number of pairs

- alpha:

  Level of significance

- Test:

  Name of the statistical test

- alternative:

  Direction of the alternative hypothesis

- Power:

  Exact power

The attribute `margin` holds the margin.

## Details

The numbers of the two kinds of discordant pairs and of concordant pairs
follow the multinomial distribution with `N` trials and probabilities
`p10`, `p01` and `1 - p10 - p01`, and the power is the probability of
the rejection region of
[`PairedRR`](https://gosukehommaex.github.io/bbssr/reference/PairedRR.md)
under this distribution, computed by complete enumeration. The
probabilities of the concordant pairs with both and with neither
treatment responding do not enter, since no test depends on how the
concordant pairs split.

## Examples

``` r
# Average exact power of the asymptotic McNemar test with 25 pairs over the alternatives
# of Berger and Sidik (2003, Section 5.1)
alt <- do.call(rbind, lapply(seq(0.025, 0.475, by = 0.05), function(q) {
  data.frame(p01 = q, p10 = seq(q + 0.05, 1 - q, by = 0.05))
}))
pw <- PairedPower(alt$p10, alt$p01, N = 25, alpha = 0.05, Test = 'McNemar')
mean(pw$Power)
#> [1] 0.6181377
```
