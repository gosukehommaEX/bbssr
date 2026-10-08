# bbssr: Blinded Sample Size Re-Estimation for Binary Endpoints

Tools for blinded sample size re-estimation (BSSR) in two-arm clinical
trials with binary endpoints, together with the exact power and sample
size calculations the re-estimation relies on. Seven tests are
supported, each available with a one-sided or a two-sided alternative.
Two of them also test non-inferiority with a margin on the scale of the
risk difference or the risk ratio, and the exact unconditional tests can
be combined with the Berger-Boos procedure.

## Main functions

- [`BinaryRR`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md):

  Rejection region of an exact test

- [`BinaryPower`](https://gosukehommaex.github.io/bbssr/reference/BinaryPower.md):

  Exact power at a given sample size

- [`BinarySampleSize`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md):

  Sample size attaining a target power

- [`BinaryPowerBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md):

  Power and final sample size of a BSSR design

- [`BinaryTypeIErrorBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md):

  Type I error rate of a BSSR design over the common response
  probability

- [`BinaryAlphaAdjBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryAlphaAdjBSSR.md):

  Adjusted significance level that controls the type I error rate of a
  BSSR design

- [`BinaryCondRejectBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryCondRejectBSSR.md):

  Rejection probabilities of a BSSR design given the pooled numbers of
  responders, and its type I error rate decomposed by the interim
  outcome

- [`BinaryGridBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryGridBSSR.md):

  Power and final sample size of several BSSR designs over common
  scenarios, as one data frame

- [`BinaryBSSR`](https://gosukehommaex.github.io/bbssr/reference/BinaryBSSR.md):

  Sample size re-estimation from observed interim data

## Reuse of computed p-values

The p-values of a test over the outcome grid depend only on the sample
sizes and the test, and a rejection region at any level follows from
them by comparison. The functions of the package therefore compute each
matrix of p-values once and keep it for the rest of the R session, which
makes a second evaluation of a related design, or a search over the
significance level, much faster than the first evaluation. The stored
matrices are discarded when they would occupy more than about 100 MB,
and `options(bbssr.cache = FALSE)` turns the reuse off. The results are
the same either way.

## See also

Useful links:

- <https://gosukehommaex.github.io/bbssr/>

- <https://github.com/gosukehommaEX/bbssr>

- Report bugs at <https://github.com/gosukehommaEX/bbssr/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
