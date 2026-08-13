# bbssr: Blinded Sample Size Re-Estimation for Binary Endpoints <img src="man/figures/bbssr_sticker.png" align="right" height="139" />

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/bbssr)](https://CRAN.R-project.org/package=bbssr)
[![R-CMD-check](https://github.com/gosukehommaEX/bbssr/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/bbssr/actions)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/bbssr)](https://cranlogs.r-pkg.org/badges/grand-total/bbssr)
[![downloads](https://cranlogs.r-pkg.org/badges/bbssr)](https://cranlogs.r-pkg.org/badges/bbssr)
<!-- badges: end -->

## Overview

A trial with a binary endpoint needs an assumed response probability in each group before
the first patient is enrolled. If the pooled response probability turns out to differ from
what was assumed, the trial ends up either underpowered or larger than it needed to be.
Blinded sample size re-estimation looks at the pooled number of responders partway through
the trial and adjusts the remaining enrolment, without ever splitting the data by treatment
group.

`bbssr` covers the whole chain of calculations that this requires, from the rejection
region of an exact test to the enrolment decision at an interim analysis.

## Version 2.0.0

Version 2 corrects the p-value of the exact unconditional tests, adds two-sided
alternatives and the Berger-Boos procedure, introduces `BinaryBSSR()` for re-estimating
from observed interim data, and moves the inner loops to C++. The `weighted` argument of
`BinaryPowerBSSR()` has been removed and `BinaryPower()` now returns a data frame. See
[NEWS.md](NEWS.md) for the full list and for what to change in code written for version 1.

## Installation

```r
install.packages("bbssr")
```

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/bbssr")
```

A C++ compiler is required to install from source.

## The five tests

| `Test` | Description | Exact level |
|---|---|---|
| `'Chisq'` | Pearson chi-squared, no continuity correction | no |
| `'Fisher'` | Fisher exact test, conditional on the margin | yes |
| `'Fisher-midP'` | Mid-p variant of the Fisher exact test | no |
| `'Z-pool'` | Exact unconditional test, Z statistic with pooled variance | yes |
| `'Boschloo'` | Exact unconditional test, Fisher p-value as ordering statistic | yes |

Only the three exact tests hold the type I error rate below the nominal level for every
value of the nuisance parameter. The Boschloo p-value never exceeds the Fisher p-value, so
the Boschloo test is uniformly the more powerful of the two.

Every test accepts `alternative = 'greater'` or `alternative = 'two.sided'`. For the
conditional tests, `tsmethod` selects between the `'minlike'` convention of
`stats::fisher.test` and the `'central'` convention that doubles the smaller tail.

## Quick start

### Power and sample size

```r
library(bbssr)

BinaryPower(p1 = 0.6, p2 = 0.3, N1 = 40, N2 = 40, alpha = 0.025, Test = 'Fisher')
#> Exact power for a two-arm trial with a binary endpoint
#> 
#>   Test         : Fisher
#>   Alternative  : greater
#>   Sample sizes : N1 = 40, N2 = 40
#>   Alpha        : 0.025
#> 
#>   p1  p2  Power
#>  0.6 0.3 0.7248

BinarySampleSize(p1 = 0.6, p2 = 0.3, r = 1, alpha = 0.025,
                 tar.power = 0.8, Test = 'Fisher')
#> Sample size for a two-arm trial with a binary endpoint
#> 
#>   Test             : Fisher
#>   Alternative      : greater
#>   Response rates   : p1 = 0.6, p2 = 0.3
#>   Allocation ratio : 1 to 1
#>   Alpha            : 0.025
#>   Target power     : 0.8
#> 
#>   Required sample size: N1 = 48, N2 = 48, total N = 96
#>   Attained power      : 0.8005
```

The exact power is not monotone in the sample size. Between 42 and 53 patients per group
the power of the test above rises from 0.7638 to 0.8540 but drops twice on the way, at 43
and at 53, because the set of attainable significance levels changes with every increment.
The search evaluates the exact power at each candidate rather than inverting a smooth
approximation.

### Rejection regions

```r
RR <- BinaryRR(N1 = 10, N2 = 10, alpha = 0.025, Test = 'Boschloo')
RR
#> Rejection region for a two-arm trial with a binary endpoint
#> 
#>   Test          : Boschloo
#>   Alternative   : greater
#>   Sample sizes  : N1 = 10, N2 = 10
#>   Alpha         : 0.025
#>   Grid points   : 100
#>   Berger-Boos   : not used
#> 
#>   Rejected outcomes: 23 of 121
#> 
#>       x2=0 x2=1 x2=2 x2=3 x2=4 x2=5 x2=6 x2=7 x2=8 x2=9 x2=10
#> x1=0  .    .    .    .    .    .    .    .    .    .    .    
#> x1=1  .    .    .    .    .    .    .    .    .    .    .    
#> x1=2  .    .    .    .    .    .    .    .    .    .    .    
#> x1=3  .    .    .    .    .    .    .    .    .    .    .    
#> x1=4  X    .    .    .    .    .    .    .    .    .    .    
#> x1=5  X    .    .    .    .    .    .    .    .    .    .    
#> x1=6  X    X    .    .    .    .    .    .    .    .    .    
#> x1=7  X    X    X    .    .    .    .    .    .    .    .    
#> x1=8  X    X    X    X    .    .    .    .    .    .    .    
#> x1=9  X    X    X    X    X    .    .    .    .    .    .    
#> x1=10 X    X    X    X    X    X    X    .    .    .    .    
#> 
#>   X marks rejection of the null hypothesis

plot(RR)
```

The whole grid of p-values is computed in one pass, so a power curve costs little more
than a single p-value. All 441 p-values of a 20 against 20 Boschloo grid take about 7 ms.

### Evaluating a re-estimation rule

```r
BinaryPowerBSSR(
  asmd.p1 = 0.45, asmd.p2 = 0.09,
  p = seq(0.19, 0.37, by = 0.03),
  Delta.A = 0.36, Delta.T = 0.36,
  N1 = 24, N2 = 24, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'
)
#>     p   p1   p2 power.BSSR power.TRAD  E.N
#>  0.19 0.37 0.01     0.9107     0.9565 38.2
#>  0.22 0.40 0.04     0.8221     0.8916 40.6
#>  0.25 0.43 0.07     0.7848     0.8442 43.5
#>  0.28 0.46 0.10     0.7727     0.8070 46.6
#>  0.31 0.49 0.13     0.7718     0.7800 49.4
#>  0.34 0.52 0.16     0.7747     0.7590 51.9
#>  0.37 0.55 0.19     0.7783     0.7388 53.9
```

The variance of a binary endpoint is largest at one half, so at a fixed risk difference the
required sample size grows as the pooled probability approaches that value. The fixed
design, sized at a pooled probability of 0.27, is over-powered below it and under-powered
above it. The re-estimation stays nearer the target in both directions, and `E.N` reports
what it costs against the 48 patients of the fixed design.

### Re-estimating from observed interim data

```r
BinaryBSSR(n1 = 12, n2 = 12, S = 8, Delta.A = 0.36, r = 1,
           alpha = 0.025, tar.power = 0.8, Test = 'Z-pool')
#> Blinded sample size re-estimation from interim data
#> 
#>   Test             : Z-pool
#>   Alternative      : greater
#>   Design rule      : unrestricted
#>   Assumed effect   : 0.36
#>   Alpha            : 0.025, target power 0.8
#> 
#> Interim data
#>   Patients         : n1 = 12, n2 = 12, total n = 24
#>   Responders       : S = 8, blinded pooled rate = 0.3333
#>   Recovered rates  : hat.p1 = 0.5133, hat.p2 = 0.1533
#> 
#> Re-estimation
#>   Required total   : N1 = 27, N2 = 27, total N = 54
#>   Still to enrol   : group 1 = 15, group 2 = 15, total = 30
#>   Final size       : N1 = 27, N2 = 27, total N = 54
#>   Power at final N : 0.8129
```

The total number of responders is the only quantity that has to leave the database, so the
analysis stays blinded and no alpha is spent.

## Functions

| Function | Purpose |
|---|---|
| `BinaryRR()` | Rejection region of an exact test over the outcome grid |
| `BinaryPower()` | Exact power at a given sample size |
| `BinarySampleSize()` | Smallest sample size attaining a target power |
| `BinaryPowerBSSR()` | Power and expected sample size of a re-estimation rule |
| `BinaryBSSR()` | Sample size re-estimation from observed interim data |

Each returns a classed object with a `print()` method, and the first four also have a
`plot()` method.

## Design rules

Under the **unrestricted** rule (`restricted = FALSE`) the re-estimated sample size is used
as it stands, so the trial may end up smaller than planned. Under the **restricted** rule
(`restricted = TRUE`) the planned sample size acts as a floor, so the trial can only grow.
The two rules give the same answer whenever the interim data call for more patients than
planned.

## Vignettes

- [`vignette("bbssr-introduction")`](https://htmlpreview.github.io/?https://github.com/gosukehommaEX/bbssr/blob/main/doc/bbssr-introduction.html)
- [`vignette("bbssr-statistical-methods")`](https://htmlpreview.github.io/?https://github.com/gosukehommaEX/bbssr/blob/main/doc/bbssr-statistical-methods.html)
- [`vignette("bbssr-interim-reestimation")`](https://htmlpreview.github.io/?https://github.com/gosukehommaEX/bbssr/blob/main/doc/bbssr-interim-reestimation.html)
- [`vignette("bbssr-validation")`](https://htmlpreview.github.io/?https://github.com/gosukehommaEX/bbssr/blob/main/doc/bbssr-validation.html)

## Validation

The conditional tests reproduce `stats::fisher.test` to machine precision. The
unconditional tests reproduce a direct evaluation of their definition to machine precision,
and agree with `Exact` and `exact2x2` up to the difference between the grids that the
packages place over the nuisance parameter. Details are in the validation vignette.

## References

Berger, R. L. and Boos, D. D. (1994). P values maximized over a confidence set for the
nuisance parameter. *Journal of the American Statistical Association*, 89, 1012-1016.

Boschloo, R. D. (1970). Raised conditional level of significance for the 2x2-table when
testing the equality of two probabilities. *Statistica Neerlandica*, 24, 1-9.

Kieser, M. (2020). *Methods and Applications of Sample Size Calculation and Recalculation
in Clinical Trials*. Springer.

Mehrotra, D. V., Chan, I. S. F. and Berger, R. L. (2003). A cautionary note on exact
unconditional inference for a difference between two independent binomial proportions.
*Biometrics*, 59, 441-450.

## License

MIT, see [LICENSE.md](LICENSE.md).

## Author

Gosuke Homma

---

This package is intended for statisticians and clinical researchers familiar with adaptive
trial designs. For a regulatory submission, consult biostatisticians and regulatory affairs
specialists about the applicable guidelines.
