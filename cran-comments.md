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
or more strata, and with random stratum and group sizes, the rates of the re-estimation
design are estimated by simulation. It also
adds functions for the rejection probabilities given the interim outcome and for the
evaluation of several designs at once, together with further options for the
re-estimation. `NEWS.md` lists every change, including the changes in behaviour.

## Test environments

* local R installation, R 4.6.1 on Windows 11
* win-builder, R-devel (2026-10-23 r90760) and R-release (4.6.1)
* GitHub Actions, R-CMD-check on Ubuntu (devel, release, oldrel-1), macOS and Windows

## R CMD check results

0 errors | 0 warnings | 0 notes on the local installation. R-CMD-check on GitHub Actions
passes on all five platforms.

win-builder (R-devel and R-release) gives 1 NOTE:

    Possibly misspelled words in DESCRIPTION:
      Blackwelder (14:20)
      Farrington (15:5)
      Friede (29:43, 30:28)
      Haenszel (20:32)
      McNemar (18:68)
      Woolson (21:17)

These are the surnames of the authors of the cited references and are spelled correctly.

TODO: re-run the checks after the additions planned for this release and update the
results above before the submission.

## Reverse dependencies

There are no reverse dependencies on CRAN.
