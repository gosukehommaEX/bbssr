## Submission

This is a major update of a package already on CRAN, from version 1.0.2 to 2.0.0.

The main reason for the release is a correction. The exact unconditional tests accumulated
the null probabilities along an arbitrary ordering of the outcomes, so outcomes sharing the
same value of the ordering statistic could receive different p-values and different
decisions. Tied values are now grouped explicitly. The release also adds two-sided
alternatives, the Berger-Boos procedure, a function for re-estimating the sample size from
observed interim data, and 'Rcpp' for the inner loops.

The user interface has changed in ways that are not backward compatible, which is why the
major version number has been raised. The `weighted` argument of `BinaryPowerBSSR()` has
been removed and `BinaryPower()` now returns a data frame rather than a numeric vector.
`NEWS.md` lists every change.

## Test environments

* local R installation, R 4.6.0 on Windows 11
* win-builder, R-devel and R-release
* macOS builder, https://mac.r-project.org/macbuilder/
* GitHub Actions, R-CMD-check on Ubuntu (devel, release, oldrel-1), macOS and Windows

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
