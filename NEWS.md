# bbssr (development version)

## New features

* `BinaryTypeIErrorBSSR()` evaluates the type I error rate of a re-estimation design and
  of the corresponding fixed-sample design over the common response probability, and
  locates the largest value. The rate is a polynomial in the response probability, so the
  largest local maxima on the grid are refined by a one-dimensional optimization.

* `BinaryAlphaAdjBSSR()` finds the adjusted nominal level that keeps the largest type I
  error rate at or below the target level (Kieser and Friede, 2000; Friede and Kieser,
  2004). The adjusted level can be applied to the final analysis only, which allows a
  bisection, or to the re-estimation as well, which uses a stepwise search.

* `BinaryCondRejectBSSR()` returns, under equal response probabilities, the rejection
  probability of a re-estimation design given the pooled numbers of responders at the
  interim analysis and in the second stage, together with the rejection probability
  given only their total. The first does not depend on the common response probability
  and the type I error rate is its binomial average, so the type I error rate is at most
  the nominal level whenever the first is, and the rate can be decomposed by the pooled
  number of interim responders to show the interim outcomes that raise it. Because the
  allocation is fixed within each stage, the first can exceed the level even for
  Fisher's exact test without re-estimation, whereas the conditional size of that test
  given the total never does. `inst/validation/kieser-2020-section-21-3.R` uses both to
  examine numerically the statement of Kieser (2020, Sect. 21.3) that the type I error
  rate is controlled for any blinded re-estimation rule when the Fisher-Boschloo test is
  applied in the analysis.

* The sample size can be re-estimated from the normal approximation instead of the exact
  power, through `ss.method = 'standard'` (formula 21.3 of Kieser, 2020) or
  `ss.method = 'null.variance'` (formula 1 of Friede and Kieser, 2004). The same options
  are available in `BinarySampleSize()` as `method`. The re-estimation can also use a
  test (`ss.Test`) or a level (`ss.alpha`) other than those of the final analysis.
  The normal approximation uses the recovered proportions before truncation to the unit
  interval, with each Bernoulli variance truncated at zero, so that the assumed effect is
  kept as in formula (2) of Friede and Kieser (2004). Truncating the proportions instead
  shrinks the assumed difference near the boundary and inflates the re-estimated size.

* `rounding` selects how an unrounded sample size becomes whole numbers: by group
  (`'group'`, the rule of earlier versions), by the rule of Friede and Kieser (2004)
  (`'friede-kieser'`), or by rounding the total (`'total'`).

* `N.min` and `N.max` bound the final total sample size of a re-estimation design, for
  example to cap it at twice the initial size or to keep patients who are enrolled but
  not yet evaluated.

* `n.interim` gives the interim sample sizes directly, as an alternative to `omega`.

* `effect = 'RR'` and `effect = 'OR'` split the blinded pooled proportion with an assumed
  risk ratio or odds ratio instead of a risk difference. Formula (21.12) of Kieser (2020)
  for the odds ratio repeats the formula for the risk ratio, so the odds ratio is
  obtained by solving the defining quadratic equation instead.

* All functions accept `alternative = 'less'`, obtained by exchanging the two groups.

* Two tests of non-inferiority are added, `Test = 'Blackwelder'` (Blackwelder, 1982),
  whose standard error uses the observed proportions, and `Test = 'Farrington-Manning'`
  (Farrington and Manning, 1990), whose standard error uses the maximum likelihood
  estimates under the null hypothesis. The argument `margin` of every function gives the
  margin on the scale of the risk difference. Without a margin the two tests are tests of
  superiority, and the Farrington-Manning test coincides with `'Chisq'`. Their rejection
  regions, power and sample sizes are computed exactly, as for the other tests.

* With a margin, the normal approximation of `method = 'standard'` is formula (4) of
  Farrington and Manning (1990), and the new `method = 'alternative.variance'` gives the
  formula of Blackwelder (1982). `rounding = 'nearest'` rounds each group to the nearest
  whole number. The blinded re-estimation of Friede, Mitchell and Mueller-Velten (2007)
  follows from `Delta.A = 0` and a margin, and `BinaryTypeIErrorBSSR()` and
  `BinaryAlphaAdjBSSR()` evaluate the type I error rate on the boundary of the null
  hypothesis, parametrized by the pooled response probability.

* `inst/reproduce/reproduce-published.R` also reproduces Tables I and II and the first
  example of Farrington and Manning (1990), Table 3 and the examples of Blackwelder
  (1982), and Tables 2 and 3 and Sections 5 and 6 of Friede, Mitchell and Mueller-Velten
  (2007). A difference from a published value counts as explained only when it is within
  the tolerance recorded with its reason.

* `ref.pvalue = TRUE` refines the maximization over the nuisance parameter of the
  Z-pooled and Boschloo tests, in every function that computes a rejection region. The
  maximum over the grid of `n.grid` points can understate the p-value, so a test can
  exceed its level: with 32 patients per group, the one-sided Z-pooled test at the level
  0.025 has a size slightly above 0.025 on the default grid and of about 0.0233 with the
  refinement. The grid is extended by points equally spaced on the arcsine square-root
  scale, and every local maximum on the extended grid is refined by a safeguarded Newton
  iteration. The default is `FALSE`, which keeps the results of earlier versions.

* `BinaryPowerBSSR()` returns the final sample size of every interim outcome as the
  attribute `reestimation` and the distribution of the final sample size as the
  attribute `N.dist`, and `summary()` reports its standard deviation and quantiles.

* `inst/reproduce/reproduce-published.R` reproduces Table I and Section 5 of Friede and
  Kieser (2004) and Example 21.1 of Kieser (2020).

## Changes in behaviour

* `BinarySampleSize()` stops with an error when the direction of `p1 - p2` contradicts
  the alternative. It previously searched without end in this case.

* `BinaryBSSR()` no longer requires `Delta.A` to lie in (0, 1), since negative risk
  differences and ratios are now allowed, and under the restricted rule it never returns
  a final size below the observed interim size.

* The p-values of each test are stored instead of its rejection regions, so a region at
  any level is obtained without recomputing them.

## Performance

* `BinaryAlphaAdjBSSR()` evaluates the type I error rate at a new level first at the grid
  point where the largest rate of the last failing level was found, and rejects the
  level at once if the rate exceeds the target level there. The largest rate over the
  grid is at least this value, so the adjusted levels are unchanged, and the evaluation
  over the whole grid and its refinement are spared for most failing levels, which
  matters most with `adjust = 'both'`.

* The exact sample size search of the re-estimation runs once for all recovered pairs of
  proportions and obtains the rejection region of each candidate sample size only once,
  so the per-call overhead of `BinarySampleSize()` is paid once per candidate size rather
  than once per pair. Each pair receives the same sample size as a separate search would
  give.

* A rejection region depends only on the sample sizes, the level and the test, yet the
  sample size search of `BinarySampleSize()` and the re-estimation of `BinaryPowerBSSR()`
  recomputed it for every candidate sample size and every interim outcome. Each region is
  now computed once and kept for the rest of the R session. Over a design grid of five
  assumed effects and nine interim fractions this replaces more than ten thousand
  evaluations of `BinaryRR()` by fewer than a hundred. The stored regions
  are discarded when they would occupy more than about 100 MB, and
  `options(bbssr.cache = FALSE)` turns the reuse off.

* `BinaryPowerBSSR()` sums the conditional power of the second stage in compiled code.
  Each column of the final rejection region is stored as runs of rejected outcomes, so
  the conditional probability over a column is a difference of binomial distribution
  functions rather than a sum over every cell. The results agree with version 2.0.0 up
  to rounding in the last digits.

# bbssr 2.0.0

This is a major release. It corrects the p-value of the exact unconditional tests, adds
two-sided alternatives, the Berger-Boos procedure and a function for re-estimating the
sample size from observed interim data, and moves the inner loops to compiled code. The
user interface has changed in ways that are not backward compatible, which is why the
major version number has been raised.

## Corrections

* The exact unconditional tests accumulated the null probabilities along an arbitrary
  ordering of the outcomes. Outcomes sharing the same value of the ordering statistic were
  therefore split, and could receive different p-values. With `N1 = N2 = 7` the outcomes
  `(x1, x2) = (5, 1)` and `(6, 2)` both have a Fisher p-value of 2/39 and cannot be
  distinguished by the Boschloo test, yet the old code assigned them 0.0216 and 0.0287 and
  rejected at the first but not the second. Tied values are now grouped explicitly and
  every member of a group receives the tail probability accumulated up to the last member
  of that group. The rejection region no longer depends on how the sorting routine breaks
  ties.

* The Z-pooled branch indexed the matrix of binomial probabilities by the responder count
  rather than by the position of that count among the retained values. The two coincided
  in the configurations reached by the old code, so the results were correct, but the
  indexing broke as soon as the retained counts were not contiguous.

* The Boschloo branch discarded outcomes whose Fisher p-value was not below the smallest
  Fisher p-value among outcomes with a non-positive risk difference. The truncation was
  sound but had no natural counterpart for a two-sided alternative, and it has been
  removed. The whole outcome grid is now scanned, which the move to compiled code makes
  affordable.

* `pmin(1, x)` with a scalar first argument drops the `dim` attribute of `x`. This affected
  the matrix of p-values returned for a two-sided chi-squared test.

* The allocation ratio was not preserved by the re-estimation when `r` was not 1. The
  interim size of group 1 was `ceiling(omega N1)` and the second-stage size was
  `ceiling(r N22)`, so the final size of group 1 was rounded up three separate times and
  did not equal `ceiling(r)` times the final size of group 2. With `r = 2`, `N2 = 33`
  and `omega = 0.8` the interim allocation was 53 to 27 rather than 54 to 27, and a
  re-estimated control group of 40 gave 79 rather than 80 in group 1. Both the interim and
  the final size of group 1 are now obtained from the size of group 2 by a single
  `ceiling(r ...)`, and the second stage follows as the difference. The allocation ratio is
  therefore exact whenever `r` is a whole number, and as close as whole numbers allow
  otherwise. Results for `r = 1` are unchanged. `BinaryBSSR()` brings group 1 to
  `ceiling(r N2.final)` for the same reason, which also corrects an imbalance that is
  already present in the observed interim data.

## New features

* All five tests accept `alternative = 'greater'` (the default, and the behaviour of
  version 1) or `alternative = 'two.sided'`. For the conditional tests the two-sided
  p-value follows the convention selected by `tsmethod`, either `'minlike'`, which matches
  `stats::fisher.test`, or `'central'`. The chi-squared and Z-pooled tests order outcomes
  by the absolute value of the Z statistic when the alternative is two-sided.

* The number of grid points used to search over the nuisance parameter of the unconditional
  tests is now the argument `n.grid`, which defaults to 100 and reproduces the grid of
  version 1. A finer grid gives larger p-values and a more accurate test.

* The Berger-Boos procedure is available through the argument `bb.gamma`. A positive value
  restricts the search over the nuisance parameter to an exact Clopper-Pearson confidence
  interval computed from the observed number of responders, and adds `bb.gamma` to the
  resulting p-value. The default of 0 disables the procedure.

* `BinaryBSSR()` re-estimates the sample size from the blinded data of a trial that is
  under way. It takes the interim sample size of each group and the pooled number of
  responders, and returns the number of patients still to enrol in each group. This
  complements `BinaryPowerBSSR()`, which evaluates a re-estimation rule at the planning
  stage.

* Every function returns a classed object with a `print()` method, and `BinaryRR()`,
  `BinaryPower()`, `BinarySampleSize()` and `BinaryPowerBSSR()` also have a `plot()`
  method. The classes are `bbssr_rr`, `bbssr_power`, `bbssr_samplesize`,
  `bbssr_powerbssr` and `bbssr_bssr`.

* The `plot()` methods accept `main`, `sub`, `xlab`, `ylab`, `base_size` and `colours`, so
  the title, the axis labels, the font size and the palette can be replaced without
  editing the returned object. A label left at `NULL` keeps the text built from the
  object and a label set to `NA` is dropped. Where a legend is drawn, `legend.title` and
  `legend.labels` rename it, and where a vertical scale is continuous, `ylim` sets its
  range through `coord_cartesian()`, which hides rather than removes the points outside
  the range. The reference lines are no longer fixed: `plot.bbssr_powerbssr()` and
  `plot.bbssr_samplesize()` take `ref.line`, and `plot.bbssr_samplesize()` also takes
  `ref.line.N2`. The methods that draw a curve also take `show.points`, which removes the
  marker at each evaluated point and leaves the curve alone, so a fine grid stays legible.
  This matters when `Delta.T` is 0, since the two power columns of
  `BinaryPowerBSSR()` are then rejection probabilities under the null hypothesis and the
  line belongs at the level of significance rather than at the target power. The defaults
  reproduce the earlier plots exactly.

* `BinaryPowerBSSR()` reports the expected total sample size in a new column `E.N`, and
  carries the interim sample sizes as the attributes `n1.interim` and `n2.interim`.

* Arguments are validated, and invalid sample sizes, significance levels, allocation
  ratios and grid sizes now produce an informative error.

## Performance

* The maximization over the nuisance parameter and the summation of the binomial mass over
  the rejection region are implemented in C++ through `Rcpp`.

* `BinaryPowerBSSR()` evaluates `BinarySampleSize()` once per distinct pair of recovered
  response probabilities and `BinaryRR()` once per distinct pair of final sample sizes,
  rather than once per interim outcome.

## Interface changes that are not backward compatible

* The `weighted` argument of `BinaryPowerBSSR()` has been removed, along with the weighted
  design it selected. Calls that supply it will fail with an unused argument error.

* The `asmd.p1` and `asmd.p2` arguments of `BinaryPowerBSSR()` have been removed. They
  entered the calculation only through the weights of the weighted design, so with that
  design gone they had no effect on the result. The planning assumption is carried by
  `Delta.A` together with the initial sample sizes `N1` and `N2`. Calls that supply them
  will fail with an unused argument error.

* `BinaryPower()` returns a one-row-per-scenario data frame of class `bbssr_power` rather
  than a bare numeric vector. Code that used the returned value directly should now read
  the `Power` column.

* `BinaryRR()` returns a logical matrix carrying the class `bbssr_rr` and the design
  settings as attributes. The matrix itself is unchanged, so indexing and `sum()` behave as
  before, but printing now goes through `print.bbssr_rr()`.

* `BinarySampleSize()` and `BinaryPowerBSSR()` gained columns in their output, and
  `BinarySampleSize()` returns integer rather than double sample sizes.

* `restricted` now defaults to `FALSE` in `BinaryPowerBSSR()` instead of having no default.

* `ggplot2` and `Rcpp` have moved into `Imports`, so a compiler is required to install the
  package from source.

## Documentation

* The vignettes have been rewritten and a fourth has been added.
  `bbssr-introduction` covers the whole package, `bbssr-statistical-methods` derives the
  tests and the re-estimation rules, `bbssr-interim-reestimation` follows a single trial
  from planning to the interim decision, and `bbssr-validation` compares the package
  against `stats`, `Exact` and `exact2x2` and reports timings.

* The test suite has been rewritten around reference implementations that follow the
  definitions of the tests directly, so the p-values are checked against a direct
  evaluation rather than against stored values.

# bbssr 1.0.2

* Fixed the title case of the `DESCRIPTION` file as requested by CRAN.

# bbssr 1.0.1

* Addressed the comments of the CRAN review.

# bbssr 1.0.0

* First release on CRAN.
