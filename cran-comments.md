## Submission

This is a major update of a package already on CRAN, from version 2.0.0 to 3.0.0.
Version 2.1.0 was not released.

The release adds three designs: tests of non-inferiority and equivalence for two
independent groups, those of Blackwelder (1982) and of Farrington and Manning (1990), with
a margin on the scale of the risk difference or the risk ratio; matched pairs, with the
McNemar tests and the score test of Nam (1997); and stratified designs, with the
Cochran-Mantel-Haenszel test, Cochran's test and the stratified exact conditional test of
Jung (2014). For each design it provides the exact power and sample size, and for blinded
sample size re-estimation the power, the type I error rate with a certified maximum over
the nuisance parameters and the adjusted significance level that controls it; with three
or more strata the rates of the re-estimation design are estimated by simulation. It also
adds functions for the rejection probabilities given the interim outcome and for the
evaluation of several designs at once, together with further options for the
re-estimation. `NEWS.md` lists every change, including the changes in behaviour.

## Test environments

* local R installation, R 4.6.1 on Windows 11
* win-builder, R-devel and R-release
* GitHub Actions, R-CMD-check on Ubuntu (devel, release, oldrel-1), macOS and Windows

## R CMD check results

0 errors | 0 warnings | 0 notes on the local installation. R-CMD-check on GitHub Actions
passes on all five platforms.

TODO: add the results of win-builder (R-devel and R-release) before the submission. Draft
for the NOTE that win-builder is expected to give, to be checked against its output:

    Possibly misspelled words in DESCRIPTION:
      al, Blackwelder, Boschloo, BSSR, et, Farrington, Friede, Haenszel, Kieser,
      McNemar, Mehrotra, Nam, Woolson

These are the surnames of the authors of the cited references, the abbreviation BSSR,
defined in the Description as blinded sample size re-estimation, and "et al.". They are
spelled correctly.

## Reverse dependencies

There are no reverse dependencies on CRAN.
