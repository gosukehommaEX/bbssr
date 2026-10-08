## Submission

This is a minor update of a package already on CRAN, from version 2.0.0 to 2.1.0.

The release adds two tests of non-inferiority, those of Blackwelder (1982) and of
Farrington and Manning (1990), with a margin on the scale of the risk difference or the
risk ratio, and the blinded sample size re-estimation of non-inferiority trials. It also
adds functions for the exact type I error rate of a re-estimation design, with a certified
maximum over the nuisance parameter, for the adjusted significance level that controls it,
for the rejection probabilities given the interim outcome and for the evaluation of several
designs at once, together with further options for the re-estimation. `NEWS.md` lists every
change, including the changes in behaviour.

## Test environments

* local R installation, R 4.6.0 on Windows 11
* win-builder, R-devel and R-release
* GitHub Actions, R-CMD-check on Ubuntu (devel, release, oldrel-1), macOS and Windows

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
