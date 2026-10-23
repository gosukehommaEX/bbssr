# Changelog

## bbssr 3.0.0

Version 2.1.0 was not released, and its changes are included in this
version.

### New features

- Stratified designs.
  [`StratRR()`](https://gosukehommaex.github.io/bbssr/reference/StratRR.md),
  [`StratPower()`](https://gosukehommaex.github.io/bbssr/reference/StratPower.md),
  [`StratSampleSize()`](https://gosukehommaex.github.io/bbssr/reference/StratSampleSize.md),
  [`StratPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/StratPowerBSSR.md),
  [`StratBSSR()`](https://gosukehommaex.github.io/bbssr/reference/StratBSSR.md),
  [`StratTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/StratTypeIErrorBSSR.md)
  and
  [`StratAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/StratAlphaAdjBSSR.md)
  provide the fixed-design and re-estimation calculations for a trial
  stratified into two or more strata with fixed group sizes in each
  stratum. Five tests are available, all based on the number of
  responders of group 1 summed over the strata: the
  Cochran-Mantel-Haenszel test (`'CMH'`) and Cochran’s test
  (`'Cochran'`), each also with the continuity correction of
  [`stats::mantelhaen.test()`](https://rdrr.io/r/stats/mantelhaen.test.html)
  (`'CMH-cc'`, `'Cochran-cc'`), and the stratified exact conditional
  test of Jung (2014) (`'Strat-Fisher'`), with the two-sided conventions
  of the Fisher exact test. The exact power is computed for any number
  of strata and any response probabilities of the strata. The normal
  approximation to the total sample size is formula (2.8) of Woolson et
  al. (1986), and formula (2.3) of Nam (1992) for the tests with a
  continuity correction, and the exact search follows
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  within the family of designs given by the fractions of the strata and
  the allocation ratio. The re-estimation divides the blinded pooled
  response probability of every stratum between the groups with an
  assumed effect common to the strata, on the scale of the risk
  difference, the risk ratio or the odds ratio. With two strata the
  power is computed exactly, the type I error rate is a polynomial in
  the response probabilities of the two strata whose largest value over
  a rectangle is certified by subdividing its bivariate Bernstein
  coefficients, and the significance level can be adjusted. With three
  or more strata the power and the type I error rate of the
  re-estimation design are estimated by simulation, with Monte Carlo
  standard errors. The sample sizes of Woolson et al. (1986, Table 1)
  and Nam (1992, Tables 1 and 2), and the p-values of the example, the
  exact sample sizes of Table 4 and the design example of Jung (2014)
  are compared with the published values in
  `inst/reproduce/reproduce-published.R`, and the p-values agree with
  [`stats::mantelhaen.test()`](https://rdrr.io/r/stats/mantelhaen.test.html).
  A new vignette,
  [`vignette('bbssr-stratified')`](https://gosukehommaex.github.io/bbssr/articles/bbssr-stratified.md),
  describes the functions.

- Paired binary data.
  [`PairedRR()`](https://gosukehommaex.github.io/bbssr/reference/PairedRR.md),
  [`PairedPower()`](https://gosukehommaex.github.io/bbssr/reference/PairedPower.md),
  [`PairedSampleSize()`](https://gosukehommaex.github.io/bbssr/reference/PairedSampleSize.md),
  [`PairedPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/PairedPowerBSSR.md),
  [`PairedBSSR()`](https://gosukehommaex.github.io/bbssr/reference/PairedBSSR.md),
  [`PairedTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/PairedTypeIErrorBSSR.md)
  and
  [`PairedAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/PairedAlphaAdjBSSR.md)
  provide the fixed-design and re-estimation calculations for a
  matched-pairs design, such as a cross-over trial, with the
  probabilities `p10` and `p01` of the two kinds of discordant pairs.
  Five tests are available: the asymptotic McNemar test (`'McNemar'`),
  the exact conditional McNemar test (`'McNemar-exact'`), its mid-p
  version (`'McNemar-midP'`), the exact unconditional test of Suissa and
  Shuster (`'McNemar-uncond'`), with the Berger-Boos procedure and the
  refinement of `ref.pvalue`, and the score test of Nam and Tango
  (`'Paired-Score'`), which also tests non-inferiority with a margin and
  equivalence. The normal approximation to the number of pairs is
  formula (9) of Nam (1997), formula (3) of Connor (1987) for a margin
  of 0, and the exact search follows
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md).
  The re-estimation uses the blinded number of discordant pairs at the
  interim analysis, with the assumed difference kept, and re-estimates
  the number of pairs by the normal approximation, by the exact search,
  or by the conditional power rule of Schepers et al. (2026), which
  takes the smallest number of pairs at which the conditional power
  given the interim number of discordant pairs reaches the target, up to
  `N.max`. Power and type I error rate are computed exactly, and the
  largest type I error rate over the probability of a discordant pair on
  the null boundary is certified as for two groups. The p-values of the
  examples of Berger and Sidik (2003), Fagerland et al. (2013),
  Nam (1997) and Tango (1998), the sample sizes of Connor (1987,
  Table 1) and Nam (1997, Table 3), and the exact sizes and average
  powers of Berger and Sidik (2003, Tables 4 and 5) are reproduced in
  the tests. A new vignette,
  [`vignette('bbssr-paired')`](https://gosukehommaex.github.io/bbssr/articles/bbssr-paired.md),
  describes the functions.

- Tests of equivalence. `alternative = 'equivalence'` in
  [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md),
  [`BinaryPower()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPower.md),
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md),
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md),
  [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md),
  [`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
  and
  [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  tests the null hypothesis `p1 - p2 <= -lower` or `p1 - p2 >= upper`
  against `-lower < p1 - p2 < upper` by two one-sided tests of
  Blackwelder or of Farrington and Manning, each at the level `alpha`.
  The margins are given as `margin = c(lower, upper)`, or as a single
  value for both. With `margin.scale = 'RR'` they are ratios with
  `0 < lower < 1 < upper`, and a single ratio `R` gives the pair of `R`
  and `1 / R`, the smaller as the lower margin. The rejection region is
  the intersection of the two one-sided regions, and the p-value of an
  outcome is the larger of its two one-sided p-values. The normal
  approximation to the sample size finds the size at which the sum of
  the approximate powers of the two one-sided tests, minus 1, reaches
  the target power, and the exact search starts from it. The type I
  error rate of a re-estimation design is evaluated on both boundaries
  of the null hypothesis, its largest value is certified on each, and
  the adjusted level controls both.
  [`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
  then stacks the two boundaries with the column `boundary`, and
  [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  reports the boundary of the largest rate in the column `boundary.adj`.

- [`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
  evaluates the type I error rate of a re-estimation design and of the
  corresponding fixed-sample design over a grid of the common response
  probability, by default `seq(0, 1, by = 0.005)`, and locates the
  largest value. The rate is a polynomial in the response probability.
  With the default `maximize = 'certified'` the polynomial is written in
  the Bernstein basis over the interval spanned by the grid, restricted
  with a margin to the part of the null boundary on which both response
  probabilities lie in the unit interval. There it does not exceed its
  largest coefficient, and the interval is halved by the algorithm of de
  Casteljau until the bound on each piece exceeds the largest value
  found by at most 1e-12. The largest value over the whole interval is
  therefore found whether or not the grid comes close to it, and the
  bound is returned in the column `bound` of the attribute `max`. If 1e5
  subdivisions do not suffice, a warning is given and the bound remains
  valid. `maximize = 'refined'` refines the three largest local maxima
  on the grid by a one-dimensional optimization, which usually but not
  always finds the maximum, and `maximize = 'grid'` takes the largest
  value on the grid.

- [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  lowers the nominal level until the largest type I error rate does not
  exceed the target level (Kieser and Friede, 2000; Friede and Kieser,
  2004). The adjusted level can be applied to the final analysis only,
  which allows a bisection, or to the re-estimation as well, which uses
  a stepwise search. With the default `maximize = 'certified'` a level
  is accepted only when the upper bound of the type I error rate over
  the interval spanned by `theta` does not exceed the target level, and
  the columns `max.TIE.bound` and `max.TIE.adj.bound` hold the bounds at
  the nominal and at the adjusted level. The reported level found by the
  bisection lies between the largest p-value rejected and the smallest
  p-value not rejected at the level where the search ends, and has at
  most six significant digits where possible. It therefore gives the
  same rejection regions whether a p-value is rejected when it is below
  the level by more than the tolerance used throughout the package,
  below the level, or at most equal to it, and it can be used as
  printed.

- [`BinaryCondRejectBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryCondRejectBSSR.md)
  returns, under equal response probabilities, the rejection probability
  of a re-estimation design given the pooled numbers of responders at
  the interim analysis and in the second stage, together with the
  rejection probability given only their total. The first does not
  depend on the common response probability and the type I error rate is
  its binomial average, so the type I error rate is at most the nominal
  level whenever the first is, and the rate can be decomposed by the
  pooled number of interim responders to show the interim outcomes that
  raise it. Because the allocation is fixed within each stage, the first
  can exceed the level even for Fisher’s exact test without
  re-estimation, whereas the conditional size of that test given the
  total never does.

- [`BinaryGridBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryGridBSSR.md)
  evaluates a set of re-estimation designs, given as the rows of a data
  frame, over common true pooled response probabilities and returns one
  tidy data frame with the power of each design and of its fixed-sample
  counterpart and the distribution of the final sample size. The initial
  sample sizes can be given or obtained from a planning proportion with
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md),
  the largest type I error rates can be added, and
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) compare the
  designs.

- The sample size can be re-estimated from the normal approximation
  instead of the exact power, through `ss.method = 'standard'` (formula
  21.3 of Kieser, 2020) or `ss.method = 'null.variance'` (formula 1 of
  Friede and Kieser, 2004). The same options are available in
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  as `method`. The re-estimation can also use a test (`ss.Test`) or a
  level (`ss.alpha`) other than those of the final analysis. The normal
  approximation uses the recovered proportions before truncation to the
  unit interval, with each Bernoulli variance truncated at zero, so that
  the assumed effect is kept as in formula (2) of Friede and Kieser
  (2004). Truncating the proportions instead shrinks the assumed
  difference near the boundary and inflates the re-estimated size.

- `search` selects how the exact sample size search of
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  and of the exact re-estimation chooses the size of group 2.
  `'crossing'` (the default, and the search of earlier versions) steps
  from the normal approximation to a size that attains the target power
  while the size one unit smaller does not. `'smallest'` returns the
  smallest size that attains the target power, and `'stable'` the
  smallest size from which every size up to a limit attains it. The
  exact power is not monotone in the sample size, so the three can
  differ: for Fisher’s exact test with response probabilities 0.6 and
  0.3, a one-sided level of 0.025 and a target power of 0.85 they give
  52, 52 and 56 patients per group. The limit is set by `search.limit`,
  by default the larger of twice the normal approximation and the normal
  approximation plus 50, and is returned as the attribute `search.limit`
  of
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  and
  [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md)
  and as the column `N2.limit` of the attribute `reestimation`.

- `rounding` selects how an unrounded sample size becomes whole numbers:
  by group (`'group'`, the rule of earlier versions), by the rule of
  Friede and Kieser (2004) (`'friede-kieser'`), by rounding the total
  (`'total'`), or by rounding each group to the nearest whole number
  (`'nearest'`).

- `N.min` and `N.max` bound the final total sample size of a
  re-estimation design, for example to cap it at twice the initial size
  or to keep patients who are enrolled but not yet evaluated.

- `n.interim` gives the interim sample sizes directly, as an alternative
  to `omega`.

- `effect = 'RR'` and `effect = 'OR'` split the blinded pooled
  proportion with an assumed risk ratio or odds ratio instead of a risk
  difference. Formula (21.12) of Kieser (2020) for the odds ratio
  repeats the formula for the risk ratio, so the odds ratio is obtained
  by solving the defining quadratic equation instead.

- All functions accept `alternative = 'less'`, obtained by exchanging
  the two groups.

- `tsmethod = 'blaker'` adds a third two-sided convention for the
  Fisher, Fisher mid-p and Boschloo tests. The tables are ordered by the
  smaller of their two one-sided tail probabilities, formula (2) of
  Mehrotra, Chan and Berger (2003). For the exact Fisher p-value this is
  the convention `'blaker'` of `exact2x2`. The exact Fisher p-value
  never exceeds that of `'central'` and equals that of `'minlike'` when
  the groups are of equal size. The Boschloo test orders the outcomes by
  the exact Fisher p-value of the selected convention. For
  `'Fisher-midP'`, the tables tied with the observed one in the ordering
  of `'minlike'` or `'blaker'`, the observed one included, contribute
  half of their probability (Fay and Hunsberger, 2021, Section 9).

- Two tests of non-inferiority are added, `Test = 'Blackwelder'`
  (Blackwelder, 1982), whose standard error uses the observed
  proportions, and `Test = 'Farrington-Manning'` (Farrington and
  Manning, 1990), whose standard error uses the maximum likelihood
  estimates under the null hypothesis. The argument `margin` of every
  function gives the margin on the scale of the risk difference. Without
  a margin the two tests are tests of superiority, and the
  Farrington-Manning test coincides with `'Chisq'`. Their rejection
  regions, power and sample sizes are computed exactly, as for the other
  tests.

- With a margin, the normal approximation of `method = 'standard'` is
  formula (4) of Farrington and Manning (1990), and the new
  `method = 'alternative.variance'` gives the formula of Blackwelder
  (1982). The blinded re-estimation of Friede, Mitchell and
  Mueller-Velten (2007) follows from `Delta.A = 0` and a margin, and
  [`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
  and
  [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  evaluate the type I error rate on the boundary of the null hypothesis,
  parametrized by the pooled response probability.

- The new argument `margin.scale` of
  [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md),
  [`BinaryPower()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPower.md),
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md),
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md),
  [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md),
  [`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
  and
  [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md),
  also available as a column of the designs of
  [`BinaryGridBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryGridBSSR.md),
  gives the margin on the scale of the risk ratio with
  `margin.scale = 'RR'`. The margin is then a ratio `R0 > 0`, and the
  null hypothesis is `p1 / p2 <= R0` for `alternative = 'greater'` and
  `p1 / p2 >= R0` for `'less'`. The two tests of non-inferiority refer
  statistic (5) of Farrington and Manning (1990), `hat.p1 - R0 hat.p2`
  divided by its standard error, to the standard normal distribution.
  The standard error of the Blackwelder test uses the observed
  proportions, which is Method 1 of the article, and that of the
  Farrington-Manning test the restricted maximum likelihood estimates of
  formula (13). `method = 'standard'` gives formula (8) of the article
  and `method = 'alternative.variance'` its Method 1. In the
  re-estimation designs the effects are risk ratios (`effect = 'RR'`),
  and the type I error rate is certified on the boundary `p1 = R0 p2`.
  The argument comes last and defaults to `'RD'`, so existing calls are
  unchanged.

- `inst/reproduce/reproduce-published.R` reproduces Table I and Section
  5 of Friede and Kieser (2004), Example 21.1, the values stated for
  Figures 21.1 and 21.2, Example 23.1 and Section 23.2 of Kieser (2020),
  Tables I and II (Methods 1 and 3) for the difference and for the
  relative risk and the two examples of Farrington and Manning (1990),
  Table 3 and the examples of Blackwelder (1982), Tables 2 and 3 and
  Sections 5 and 6 of Friede, Mitchell and Mueller-Velten (2007), the
  table of power and the example of Boschloo (1970), the example and
  Tables 1 and 3 of Mehrotra, Chan and Berger (2003), Example 2 of
  Berger and Boos (1994), and the example and Table 1 of Section 8 of
  Fay and Hunsberger (2021). A difference from a published value counts
  as explained only when it is within the tolerance recorded with its
  reason.

- `ref.pvalue = TRUE` refines the maximization over the nuisance
  parameter of the Z-pooled and Boschloo tests, in every function that
  computes a rejection region. The maximum over the grid of `n.grid`
  points can understate the p-value, so a test can exceed its level:
  with 32 patients per group, the one-sided Z-pooled test at the level
  0.025 has a size slightly above 0.025 on the default grid and of about
  0.0233 with the refinement. The grid is extended by points equally
  spaced on the arcsine square-root scale, and every local maximum on
  the extended grid is refined by a safeguarded Newton iteration. The
  default is `FALSE`, which keeps the results of earlier versions.

- [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  returns the final sample size of every interim outcome as the
  attribute `reestimation` and the distribution of the final sample size
  as the attribute `N.dist`, and
  [`summary()`](https://rdrr.io/r/base/summary.html) reports its
  standard deviation and quantiles.

### Changes in behaviour

- [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  stops with an error when the direction of `p1 - p2` contradicts the
  alternative. It previously searched without end in this case.

- [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md)
  no longer requires `Delta.A` to lie in (0, 1), since negative risk
  differences and ratios are now allowed, and under the restricted rule
  it never returns a final size below the observed interim size.

- The p-values of each test are stored instead of its rejection regions,
  so a region at any level is obtained without recomputing them.

- [`print()`](https://rdrr.io/r/base/print.html) of the result of
  [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  stopped with an error when the attributes of the result had been
  removed by selecting some of its rows or columns; such an object is
  now printed as a data frame. The print methods of the results for
  matched pairs and stratified designs do the same.

- Sizes obtained by rounding up a product are rounded up with a
  tolerance for floating-point error: the interim sizes
  `ceiling(omega * N2)` and `ceiling(r * ceiling(omega * N2))` of the
  re-estimation functions and the size `ceiling(r * N2)` of group 1. For
  example, `omega = 0.55` with `N2 = 100` now gives 55 interim patients
  in group 2 instead of 56, and `r = 1.1` with 50 patients in group 2
  gives 55 patients in group 1 instead of 56.

- `N.min` and `N.max` must be finite numbers. `N.max = Inf` made
  [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md)
  and
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  run without end with the default rounding.

### Performance

- [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  evaluates the type I error rate at a new level first at the common
  response probability where the largest rate of the last failing level
  was found, and rejects the level at once if the rate exceeds the
  target level there. The largest rate is at least this value, so such a
  level also fails on the largest rate, and the evaluation over the
  whole grid and the search for the maximum are spared for most failing
  levels, which matters most with `adjust = 'both'`. With
  `maximize = 'grid'` the adjusted levels are unchanged. With the other
  methods this value can lie between the grid points, so a level can
  fail that the grid would accept, which can only lower the adjusted
  level. With `maximize = 'certified'` the bisection of
  `adjust = 'test'` certifies only the level it finds, and searches
  again below that level, certifying every level that passes, if the
  certification fails.

- The exact sample size search of the re-estimation runs once for all
  recovered pairs of proportions and obtains the rejection region of
  each candidate sample size only once, so the per-call overhead of
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  is paid once per candidate size rather than once per pair. Each pair
  receives the same sample size as a separate search would give.

- A rejection region depends only on the sample sizes, the level and the
  test, yet the sample size search of
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  and the re-estimation of
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  recomputed it for every candidate sample size and every interim
  outcome. The p-values behind each region are now computed once and
  kept for the rest of the R session, and a region at any level follows
  from them. Over a design grid of five assumed effects and nine interim
  fractions this replaces more than ten thousand evaluations of
  [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)
  by fewer than a hundred. The stored p-values are discarded when they
  would occupy more than about 100 MB, and
  `options(bbssr.cache = FALSE)` turns the reuse off.

- [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  sums the conditional power of the second stage in compiled code. Each
  column of the final rejection region is stored as runs of rejected
  outcomes, so the conditional probability over a column is a difference
  of binomial distribution functions rather than a sum over every cell.
  The results agree with version 2.0.0 up to rounding in the last
  digits.

### Documentation

- The vignettes are reorganized into nine. `bbssr-introduction`,
  `bbssr-statistical-methods`, `bbssr-interim-reestimation` and
  `bbssr-validation` cover the new tests and functions, and five are
  added: `bbssr-reestimation-rules` on the options of the re-estimation,
  `bbssr-type1-error` on the type I error rate and the adjusted level,
  `bbssr-non-inferiority` on the tests with a margin,
  `bbssr-design-grid` on the comparison of designs with
  [`BinaryGridBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryGridBSSR.md),
  and `bbssr-published-figures`, which redraws fifteen figures of Friede
  and Kieser (2004), Kieser (2020), Friede, Mitchell and
  Mueller-Velten (2007) and Boschloo (1970) from values computed with
  the package. Those values are computed by
  `inst/reproduce/reproduce-figures.R` and stored in
  `inst/extdata/published-figures/`, and the figures are drawn by
  `inst/reproduce/plot-figures.R`. The numbers quoted in the text are
  computed when the vignette is built.

- The help page of
  [`BinaryAlphaAdjBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md)
  and the `bbssr-type1-error` vignette describe how the type I error
  rate is controlled within a confidence interval for the nuisance
  parameter (Kieser, 2020, Section 22.1.3), with Example 23.1 of
  Kieser (2020) as the illustration.

- The `bbssr-validation` vignette checks the mid-p values of every
  two-sided convention against their definition, and the one-sided and
  `central` mid-p values against `exact2x2`.

- `README.md` is generated from `README.Rmd`, and a ‘pkgdown’ site is
  built from the help pages and the vignettes.

- `DESCRIPTION`, the help page of
  [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)
  and the vignettes cite Boschloo (1970), where the Boschloo test was
  introduced as Fisher’s test at a raised conditional level, and the
  `bbssr-statistical-methods` vignette explains why this is the same
  test as the one defined in the package.

- The help page of
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  described the exact search as returning the smallest sample size that
  attains the target power. It now describes what the search returns:
  starting from the normal approximation, a size of group 2 whose exact
  power attains the target while the size one unit smaller does not.
  Since the exact power is not monotone in the sample size, this need
  not be the smallest such size.

## bbssr 2.0.0

CRAN release: 2026-08-20

This is a major release. It corrects the p-value of the exact
unconditional tests, adds two-sided alternatives, the Berger-Boos
procedure and a function for re-estimating the sample size from observed
interim data, and moves the inner loops to compiled code. The user
interface has changed in ways that are not backward compatible, which is
why the major version number has been raised.

### Corrections

- The exact unconditional tests accumulated the null probabilities along
  an arbitrary ordering of the outcomes. Outcomes sharing the same value
  of the ordering statistic were therefore split, and could receive
  different p-values. With `N1 = N2 = 7` the outcomes
  `(x1, x2) = (5, 1)` and `(6, 2)` both have a Fisher p-value of 2/39
  and cannot be distinguished by the Boschloo test, yet the old code
  assigned them 0.0216 and 0.0287 and rejected at the first but not the
  second. Tied values are now grouped explicitly and every member of a
  group receives the tail probability accumulated up to the last member
  of that group. The rejection region no longer depends on how the
  sorting routine breaks ties.

- The Z-pooled branch indexed the matrix of binomial probabilities by
  the responder count rather than by the position of that count among
  the retained values. The two coincided in the configurations reached
  by the old code, so the results were correct, but the indexing broke
  as soon as the retained counts were not contiguous.

- The Boschloo branch discarded outcomes whose Fisher p-value was not
  below the smallest Fisher p-value among outcomes with a non-positive
  risk difference. The truncation was sound but had no natural
  counterpart for a two-sided alternative, and it has been removed. The
  whole outcome grid is now scanned, which the move to compiled code
  makes affordable.

- `pmin(1, x)` with a scalar first argument drops the `dim` attribute of
  `x`. This affected the matrix of p-values returned for a two-sided
  chi-squared test.

- The allocation ratio was not preserved by the re-estimation when `r`
  was not 1. The interim size of group 1 was `ceiling(omega N1)` and the
  second-stage size was `ceiling(r N22)`, so the final size of group 1
  was rounded up three separate times and did not equal `ceiling(r)`
  times the final size of group 2. With `r = 2`, `N2 = 33` and
  `omega = 0.8` the interim allocation was 53 to 27 rather than 54 to
  27, and a re-estimated control group of 40 gave 79 rather than 80 in
  group 1. Both the interim and the final size of group 1 are now
  obtained from the size of group 2 by a single `ceiling(r ...)`, and
  the second stage follows as the difference. The allocation ratio is
  therefore exact whenever `r` is a whole number, and as close as whole
  numbers allow otherwise. Results for `r = 1` are unchanged.
  [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md)
  brings group 1 to `ceiling(r N2.final)` for the same reason, which
  also corrects an imbalance that is already present in the observed
  interim data.

### New features

- All five tests accept `alternative = 'greater'` (the default, and the
  behaviour of version 1) or `alternative = 'two.sided'`. For the
  conditional tests the two-sided p-value follows the convention
  selected by `tsmethod`, either `'minlike'`, which matches
  [`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html), or
  `'central'`. The chi-squared and Z-pooled tests order outcomes by the
  absolute value of the Z statistic when the alternative is two-sided.

- The number of grid points used to search over the nuisance parameter
  of the unconditional tests is now the argument `n.grid`, which
  defaults to 100 and reproduces the grid of version 1. A finer grid
  gives larger p-values and a more accurate test.

- The Berger-Boos procedure is available through the argument
  `bb.gamma`. A positive value restricts the search over the nuisance
  parameter to an exact Clopper-Pearson confidence interval computed
  from the observed number of responders, and adds `bb.gamma` to the
  resulting p-value. The default of 0 disables the procedure.

- [`BinaryBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md)
  re-estimates the sample size from the blinded data of a trial that is
  under way. It takes the interim sample size of each group and the
  pooled number of responders, and returns the number of patients still
  to enrol in each group. This complements
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md),
  which evaluates a re-estimation rule at the planning stage.

- Every function returns a classed object with a
  [`print()`](https://rdrr.io/r/base/print.html) method, and
  [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md),
  [`BinaryPower()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPower.md),
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  and
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  also have a [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method. The classes are `bbssr_rr`, `bbssr_power`, `bbssr_samplesize`,
  `bbssr_powerbssr` and `bbssr_bssr`.

- The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  accept `main`, `sub`, `xlab`, `ylab`, `base_size` and `colours`, so
  the title, the axis labels, the font size and the palette can be
  replaced without editing the returned object. A label left at `NULL`
  keeps the text built from the object and a label set to `NA` is
  dropped. Where a legend is drawn, `legend.title` and `legend.labels`
  rename it, and where a vertical scale is continuous, `ylim` sets its
  range through
  [`coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html),
  which hides rather than removes the points outside the range. The
  reference lines are no longer fixed:
  [`plot.bbssr_powerbssr()`](https://gosukehommaex.github.io/bbssr/reference/plot.bbssr_powerbssr.md)
  and
  [`plot.bbssr_samplesize()`](https://gosukehommaex.github.io/bbssr/reference/plot.bbssr_samplesize.md)
  take `ref.line`, and
  [`plot.bbssr_samplesize()`](https://gosukehommaex.github.io/bbssr/reference/plot.bbssr_samplesize.md)
  also takes `ref.line.N2`. The methods that draw a curve also take
  `show.points`, which removes the marker at each evaluated point and
  leaves the curve alone, so a fine grid stays legible. This matters
  when `Delta.T` is 0, since the two power columns of
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  are then rejection probabilities under the null hypothesis and the
  line belongs at the level of significance rather than at the target
  power. The defaults reproduce the earlier plots exactly.

- [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  reports the expected total sample size in a new column `E.N`, and
  carries the interim sample sizes as the attributes `n1.interim` and
  `n2.interim`.

- Arguments are validated, and invalid sample sizes, significance
  levels, allocation ratios and grid sizes now produce an informative
  error.

### Performance

- The maximization over the nuisance parameter and the summation of the
  binomial mass over the rejection region are implemented in C++ through
  `Rcpp`.

- [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  evaluates
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  once per distinct pair of recovered response probabilities and
  [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)
  once per distinct pair of final sample sizes, rather than once per
  interim outcome.

### Interface changes that are not backward compatible

- The `weighted` argument of
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  has been removed, along with the weighted design it selected. Calls
  that supply it will fail with an unused argument error.

- The `asmd.p1` and `asmd.p2` arguments of
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  have been removed. They entered the calculation only through the
  weights of the weighted design, so with that design gone they had no
  effect on the result. The planning assumption is carried by `Delta.A`
  together with the initial sample sizes `N1` and `N2`. Calls that
  supply them will fail with an unused argument error.

- [`BinaryPower()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPower.md)
  returns a one-row-per-scenario data frame of class `bbssr_power`
  rather than a bare numeric vector. Code that used the returned value
  directly should now read the `Power` column.

- [`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)
  returns a logical matrix carrying the class `bbssr_rr` and the design
  settings as attributes. The matrix itself is unchanged, so indexing
  and [`sum()`](https://rdrr.io/r/base/sum.html) behave as before, but
  printing now goes through
  [`print.bbssr_rr()`](https://gosukehommaex.github.io/bbssr/reference/print.bbssr_rr.md).

- [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  and
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  gained columns in their output, and
  [`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
  returns integer rather than double sample sizes.

- `restricted` now defaults to `FALSE` in
  [`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
  instead of having no default.

- `ggplot2` and `Rcpp` have moved into `Imports`, so a compiler is
  required to install the package from source.

### Documentation

- The vignettes have been rewritten and a fourth has been added.
  `bbssr-introduction` covers the whole package,
  `bbssr-statistical-methods` derives the tests and the re-estimation
  rules, `bbssr-interim-reestimation` follows a single trial from
  planning to the interim decision, and `bbssr-validation` compares the
  package against `stats`, `Exact` and `exact2x2` and reports timings.

- The test suite has been rewritten around reference implementations
  that follow the definitions of the tests directly, so the p-values are
  checked against a direct evaluation rather than against stored values.

## bbssr 1.0.2

CRAN release: 2025-06-18

- Fixed the title case of the `DESCRIPTION` file as requested by CRAN.

## bbssr 1.0.1

- Addressed the comments of the CRAN review.

## bbssr 1.0.0

- First release on CRAN.
