#' Rejection Region for Two-Arm Trials with Binary Endpoints
#'
#' Provides a rejection region (RR) for two-arm trials with binary endpoints. Seven tests
#' are supported. Each can be applied with a one-sided or a two-sided alternative, and two
#' of them also test non-inferiority with a margin on the scale of the risk difference.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, \code{'Boschloo'}, \code{'Blackwelder'} or
#'   \code{'Farrington-Manning'}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default) for the one-sided alternative that the response probability of group 1
#'   exceeds that of group 2, \code{'less'} for the one-sided alternative that it falls
#'   below that of group 2, or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests, see Details. Options: \code{'minlike'} (default), \code{'central'} or
#'   \code{'blaker'}. The Boschloo test orders the outcomes by the two-sided Fisher
#'   p-value of the selected convention. Ignored for a one-sided alternative, and ignored
#'   by \code{'Chisq'}, \code{'Z-pool'}, \code{'Blackwelder'} and
#'   \code{'Farrington-Manning'}, whose two-sided versions are based on the absolute value
#'   of the Z statistic
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100. Ignored by the conditional tests
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure for the
#'   unconditional tests. The default of 0 disables the procedure. A positive value
#'   restricts the search over the nuisance parameter to an exact
#'   \code{100 (1 - bb.gamma)} percent confidence interval and adds \code{bb.gamma} to the
#'   resulting p-value. A common choice is 0.0001
#' @param margin Non-inferiority margin on the scale of the risk difference. The
#'   default of 0 gives a test of superiority. A value other than 0 tests the null
#'   hypothesis \code{p1 - p2 <= -margin} against \code{p1 - p2 > -margin} when
#'   \code{alternative} is \code{'greater'}, and \code{p1 - p2 >= margin} against
#'   \code{p1 - p2 < margin} when it is \code{'less'}. It requires
#'   \code{Test = 'Blackwelder'} or \code{'Farrington-Manning'}, see
#'   \code{\link{BinaryRR}}. A negative value tests for superiority by more than its
#'   absolute value
#' @param ref.pvalue Logical. If \code{TRUE}, the maximization over the nuisance parameter
#'   of the unconditional tests is refined between the grid points, see Details. Default
#'   is \code{FALSE}. Ignored by the conditional tests
#'
#' @return An object of class \code{bbssr_rr}, which is a logical matrix of dimension
#'   \code{(N1 + 1)} by \code{(N2 + 1)} whose entry \code{[i + 1, j + 1]} is \code{TRUE}
#'   when the null hypothesis is rejected at \code{i} responders in group 1 and \code{j}
#'   responders in group 2. The design settings are stored as attributes
#'
#' @details
#' The function supports the following seven tests:
#' \itemize{
#'   \item The Pearson chi-squared test (Chisq)
#'   \item The Fisher exact test (Fisher)
#'   \item The Fisher mid-p test (Fisher-midP)
#'   \item The Z-pooled exact unconditional test (Z-pool)
#'   \item The exact unconditional test of Boschloo (1970), which orders the outcomes by
#'     the Fisher p-value (Boschloo)
#'   \item The test of Blackwelder (1982) with the unpooled standard error (Blackwelder)
#'   \item The test of Farrington and Manning (1990) with the standard error at the
#'     restricted maximum likelihood estimates (Farrington-Manning)
#' }
#'
#' For the two-sided versions of the conditional tests, \code{'minlike'} sums the null
#' probabilities of all tables that are no more likely than the observed table, which is the
#' convention of \code{stats::fisher.test}, and \code{'central'} doubles the smaller of
#' the two one-sided tail probabilities. \code{'blaker'} orders the tables by the smaller
#' of their two one-sided tail probabilities and sums the null probabilities of all tables
#' at which this is no larger than at the observed table, formula (2) of Mehrotra, Chan
#' and Berger (2003), which is the convention \code{'blaker'} of the \pkg{exact2x2}
#' package.
#' Its Fisher p-value never exceeds that of \code{'central'}, and it equals that of
#' \code{'minlike'} when the two groups are of equal size, since the conditional
#' distribution is then symmetric. Under the mid-p correction the tables tied with the
#' observed table in the ordering of \code{'minlike'} or \code{'blaker'}, the observed
#' table included, contribute half of their probability, following the definition of the
#' mid-p value in Fay and Hunsberger (2021, Section 9). The two-sided versions of
#' \code{'Chisq'} and \code{'Z-pool'} order the outcomes by the absolute value of the Z
#' statistic.
#'
#' The unconditional tests maximize the null tail probability of an ordering statistic over
#' the common response probability, which is a nuisance parameter. Outcomes sharing the same
#' value of the ordering statistic receive the same p-value. With a positive
#' \code{bb.gamma} the maximization is restricted to a confidence interval for the
#' nuisance parameter and \code{bb.gamma} is added, following Berger and Boos (1994).
#'
#' The tail probability of each outcome is a polynomial in the common response
#' probability, and its maximum over the grid of \code{n.grid} points is a lower bound of
#' the maximum over the unit interval. A p-value computed on the grid can therefore fall
#' below the exact p-value, and the test can then exceed its nominal level. With
#' \code{ref.pvalue = TRUE} the grid is extended by points equally spaced on the arcsine
#' square-root scale, which resolve the narrow local maxima near 0 and 1, and every local
#' maximum on the extended grid is refined by a safeguarded Newton iteration between its
#' two neighbouring grid points. The grid maximum is kept when it is larger, so a refined
#' p-value is never smaller than the p-value on the grid. In the designs examined during
#' development the refined p-values agreed with a certified maximum over the unit interval
#' to within 1e-12, but the refinement is a local search and does not guarantee the
#' maximum. It takes longer than the grid alone, by roughly an order of magnitude with 500
#' patients per group.
#'
#' The alternative \code{'less'} is handled by exchanging the two groups, testing the
#' alternative \code{'greater'} and exchanging them back. Every test considered treats
#' the groups symmetrically apart from the direction of the alternative, so this gives
#' the lower-tail version of each test.
#'
#' The tests of Blackwelder and of Farrington and Manning refer the statistic
#' \code{z = (hat.p1 - hat.p2 + margin) / SE} to the standard normal distribution, which
#' tests the null hypothesis \code{p1 - p2 <= -margin} for \code{alternative = 'greater'}.
#' The standard error of the Blackwelder test uses the observed proportions. That of the
#' Farrington-Manning test uses the maximum likelihood estimates under
#' \code{p1 - p2 = -margin}, formula (12) of Farrington and Manning (1990). With
#' \code{margin = 0} the Farrington-Manning test coincides with \code{'Chisq'} and the
#' Blackwelder test is the Wald test. A standard error of 0, which arises only when both
#' observed proportions are 0 or 1, gives \code{z = Inf} or \code{-Inf} according to the
#' sign of the numerator, and \code{z = 0} when the numerator is 0 as well. These tests are
#' asymptotic, and their rejection regions and power are computed exactly, as for the
#' other tests. A margin other than 0 is available only with these two tests and a
#' one-sided alternative.
#'
#' The p-values are kept for the rest of the session and reused by later calls with the
#' same sample sizes and test, see \code{\link{bbssr-package}}.
#'
#' @references
#' Berger RL, Boos DD (1994). P values maximized over a confidence set for the nuisance
#' parameter. \emph{Journal of the American Statistical Association}, 89(427),
#' 1012-1016.
#'
#' Blackwelder WC (1982). "Proving the null hypothesis" in clinical trials.
#' \emph{Controlled Clinical Trials}, 3(4), 345-353.
#'
#' Boschloo RD (1970). Raised conditional level of significance for the 2 x 2-table when
#' testing the equality of two probabilities. \emph{Statistica Neerlandica}, 24(1), 1-9.
#'
#' Farrington CP, Manning G (1990). Test statistics and sample size formulae for
#' comparative binomial trials with null hypothesis of non-zero risk difference or
#' non-unity relative risk. \emph{Statistics in Medicine}, 9(12), 1447-1454.
#'
#' Fay MP, Hunsberger SA (2021). Practical valid inferences for the two-sample binomial
#' problem. \emph{Statistics Surveys}, 15, 72-110.
#'
#' Mehrotra DV, Chan ISF, Berger RL (2003). A cautionary note on exact unconditional
#' inference for a difference between two independent binomial proportions.
#' \emph{Biometrics}, 59(2), 441-450.
#'
#' @examples
#' # Simple example with small sample sizes
#' RR <- BinaryRR(N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
#' print(RR)
#'
#' \donttest{
#' # Two-sided Boschloo test with the Berger-Boos procedure
#' RR <- BinaryRR(N1 = 20, N2 = 10, alpha = 0.05, Test = 'Boschloo',
#'                alternative = 'two.sided', bb.gamma = 0.0001)
#' print(RR)
#' plot(RR)
#'
#' # The grid maximum can understate the p-value. Refinement removes two outcomes from the
#' # rejection region of the Z-pooled test with 32 patients per group
#' RR.grid <- BinaryRR(N1 = 32, N2 = 32, alpha = 0.025, Test = 'Z-pool')
#' RR.ref <- BinaryRR(N1 = 32, N2 = 32, alpha = 0.025, Test = 'Z-pool', ref.pvalue = TRUE)
#' sum(RR.grid) - sum(RR.ref)
#'
#' # Non-inferiority with a margin of 0.1 by the test of Farrington and Manning
#' RR <- BinaryRR(N1 = 60, N2 = 60, alpha = 0.025, Test = 'Farrington-Manning',
#'                margin = 0.1)
#' print(RR)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
BinaryRR <- function(N1, N2, alpha, Test,
                     alternative = c('greater', 'less', 'two.sided'),
                     tsmethod = c('minlike', 'central', 'blaker'),
                     n.grid = 100, bb.gamma = 0, margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  a <- check_rr_args(N1, N2, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                     margin)
  Test <- a$Test
  N1 <- a$N1
  N2 <- a$N2
  n.grid <- a$n.grid
  ref.pvalue <- a$ref.pvalue
  margin <- a$margin
  p.val <- get_pvalue(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma, ref.pvalue,
                      margin)
  RR <- (p.val %<<% alpha)
  dim(RR) <- c(N1 + 1L, N2 + 1L)
  dimnames(RR) <- list(x1 = as.character(0:N1), x2 = as.character(0:N2))
  attr(RR, 'N1') <- N1
  attr(RR, 'N2') <- N2
  attr(RR, 'alpha') <- alpha
  attr(RR, 'Test') <- Test
  attr(RR, 'alternative') <- alternative
  attr(RR, 'tsmethod') <- tsmethod
  attr(RR, 'n.grid') <- n.grid
  attr(RR, 'bb.gamma') <- bb.gamma
  attr(RR, 'margin') <- margin
  attr(RR, 'ref.pvalue') <- ref.pvalue
  attr(RR, 'p.value') <- p.val
  class(RR) <- c('bbssr_rr', 'matrix', 'array')
  RR
}
