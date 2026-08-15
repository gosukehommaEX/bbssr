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
  `ceiling(r N22)`, so the final size of group 1 was the sum of three separate roundings
  and did not equal `ceiling(r)` times the final size of group 2. With `r = 2`, `N2 = 33`
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
