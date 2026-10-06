#' Sample Size Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates the required sample size for two-arm trials with binary endpoints. Seven
#' tests are supported (see \code{\link{BinaryRR}}), each of which can be applied with a
#' one-sided or a two-sided alternative. The sample size can be obtained from the exact
#' power of the selected test or from the normal approximation.
#'
#' @param p1 True probability of responders for group 1
#' @param p2 True probability of responders for group 2
#' @param r Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2 = r:1,
#'   r > 0)
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, \code{'Boschloo'}, \code{'Blackwelder'} or
#'   \code{'Farrington-Manning'}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default), which requires \code{p1 - p2 > -margin}, \code{'less'}, which requires
#'   \code{p1 - p2 < margin}, or \code{'two.sided'}, which requires \code{p1 != p2}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests. Options: \code{'minlike'} (default) or \code{'central'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#' @param method How the sample size is obtained. \code{'exact'} (default) steps from the
#'   normal approximation, one patient at a time, to the sample size at which the exact
#'   power of \code{Test} first attains the target, see Details.
#'   \code{'standard'} uses the normal approximation with the variance under the null
#'   hypothesis for the significance term and the variance under the alternative for the
#'   power term, formula (21.3) of Kieser (2020). \code{'null.variance'} uses the variance
#'   under the null hypothesis for both terms, formula (1) of Friede and Kieser (2004).
#'   \code{'alternative.variance'} uses the variance under the alternative for both terms,
#'   as in Blackwelder (1982). See Details for a non-zero \code{margin}
#' @param rounding How an unrounded sample size from the normal approximation is turned
#'   into whole numbers. \code{'group'} (default) rounds up the size of group 2 and gives
#'   group 1 \code{ceiling(r N2)} patients. \code{'friede-kieser'} rounds up the two
#'   group sizes separately, as in Friede and Kieser (2004). \code{'total'} rounds up the
#'   total and gives group 2 \code{floor(N / (1 + r))} patients. \code{'nearest'} rounds
#'   the two group sizes separately to the nearest whole number, as in Farrington and
#'   Manning (1990). Only \code{'group'} is available with \code{method = 'exact'}
#' @param margin Non-inferiority margin on the scale of the risk difference. The
#'   default of 0 gives a test of superiority. A value other than 0 tests the null
#'   hypothesis \code{p1 - p2 <= -margin} against \code{p1 - p2 > -margin} when
#'   \code{alternative} is \code{'greater'}, and \code{p1 - p2 >= margin} against
#'   \code{p1 - p2 < margin} when it is \code{'less'}. It requires
#'   \code{Test = 'Blackwelder'} or \code{'Farrington-Manning'}, see
#'   \code{\link{BinaryRR}}. A negative value tests for superiority by more than its
#'   absolute value
#' @param ref.pvalue Logical. If \code{TRUE}, the maximization over the nuisance parameter
#'   of the unconditional tests is refined between the grid points, see
#'   \code{\link{BinaryRR}}. Default is \code{FALSE}
#'
#' @return An object of class \code{bbssr_samplesize}, a data frame with one row
#'   containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{r}{Allocation ratio to group 1}
#'   \item{alpha}{Level of significance}
#'   \item{tar.power}{Target power}
#'   \item{Test}{Name of the statistical test}
#'   \item{alternative}{Direction of the alternative hypothesis}
#'   \item{Power}{Exact power of \code{Test} at the selected sample size}
#'   \item{N1}{Required sample size of group 1}
#'   \item{N2}{Required sample size of group 2}
#'   \item{N}{Total required sample size}
#' }
#'
#' @details
#' The calculation uses a three-step approach:
#' \enumerate{
#'   \item Calculate an initial sample size from the normal approximation to the
#'     chi-squared test
#'   \item Evaluate the exact power at the initial sample size
#'   \item If the exact power attains the target power, lower the size of group 2 one
#'     unit at a time as long as the power still attains it; otherwise raise it one
#'     unit at a time until the power attains it
#' }
#'
#' The normal approximation of the first step uses \code{alpha} for a one-sided
#' alternative and \code{alpha / 2} for a two-sided alternative. Only the starting value of
#' the search is affected, so the returned sample size is exact in either case.
#'
#' The search stops where the exact power crosses the target: the returned size of
#' group 2 attains the target power and the size one unit smaller does not, unless the
#' returned size is 1. The exact power is not monotone in the sample size, so a smaller
#' size further from the normal approximation may also attain the target power, and a
#' larger size may fall short of it.
#'
#' Under \code{method = 'standard'}, \code{'null.variance'} or
#' \code{'alternative.variance'} the steps above are replaced by the closed-form normal
#' approximation and the rounding rule selected by \code{rounding}. The \code{Power}
#' column still reports the exact power of \code{Test} at the resulting sample size.
#'
#' With a non-zero \code{margin}, the difference \code{p1 - p2} in the normal
#' approximation is replaced by its distance \code{p1 - p2 + margin} from the boundary of
#' the null hypothesis, and the variance under the null hypothesis is evaluated at the
#' large sample values of the restricted maximum likelihood estimates of Farrington and
#' Manning (1990). The method \code{'standard'} then gives formula (4) of Farrington and
#' Manning (1990), which is formula (2) of Friede et al. (2007), and
#' \code{'alternative.variance'} gives the formula of Blackwelder (1982), which is formula
#' (1) of Friede et al. (2007). The exact search starts from the formula of Farrington and
#' Manning.
#'
#' @references
#' Friede T, Kieser M (2004). Sample size recalculation for binary data in internal pilot
#' study designs. \emph{Pharmaceutical Statistics}, 3(4), 269-279.
#'
#' Kieser M (2020). \emph{Methods and Applications of Sample Size Calculation and
#' Recalculation in Clinical Trials}. Springer, Cham.
#'
#' Blackwelder WC (1982). "Proving the null hypothesis" in clinical trials.
#' \emph{Controlled Clinical Trials}, 3(4), 345-353.
#'
#' Farrington CP, Manning G (1990). Test statistics and sample size formulae for
#' comparative binomial trials with null hypothesis of non-zero risk difference or
#' non-unity relative risk. \emph{Statistics in Medicine}, 9(12), 1447-1454.
#'
#' Friede T, Mitchell C, Mueller-Velten G (2007). Blinded sample size reestimation in
#' non-inferiority trials with binary endpoints. \emph{Biometrical Journal}, 49(6),
#' 903-916.
#'
#' @examples
#' # One-sided chi-squared test
#' BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
#'                  tar.power = 0.8, Test = 'Chisq')
#'
#' \donttest{
#' # Two-sided Fisher exact test
#' BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 2, alpha = 0.05,
#'                  tar.power = 0.9, Test = 'Fisher', alternative = 'two.sided')
#'
#' # Normal approximation with the variance under the null hypothesis, two-sided test
#' # at level 0.05, as in Friede and Kieser (2004)
#' BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.05, tar.power = 0.8,
#'                  Test = 'Chisq', alternative = 'two.sided',
#'                  method = 'null.variance')
#'
#' # Non-inferiority with a margin of 0.1, formula (2) of Friede et al. (2007)
#' BinarySampleSize(p1 = 0.7, p2 = 0.7, r = 1, alpha = 0.025, tar.power = 0.8,
#'                  Test = 'Farrington-Manning', method = 'standard',
#'                  rounding = 'friede-kieser', margin = 0.1)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @importFrom stats dbinom
BinarySampleSize <- function(p1, p2, r, alpha, tar.power, Test,
                             alternative = c('greater', 'less', 'two.sided'),
                             tsmethod = c('minlike', 'central'),
                             n.grid = 100, bb.gamma = 0,
                             method = c('exact', 'standard', 'null.variance',
                                        'alternative.variance'),
                             rounding = c('group', 'friede-kieser', 'total', 'nearest'),
                             margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  method <- match.arg(method)
  rounding <- match.arg(rounding)
  if (length(p1) != 1 || length(p2) != 1) stop('p1 and p2 must each be a single value')
  if (length(r) != 1 || is.na(r) || r <= 0) stop('r must be a single positive value')
  if (length(tar.power) != 1 || tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  # Validates the test, the level, the margin and the remaining arguments of the
  # rejection region
  a <- check_rr_args(1, 1, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                     margin)
  Test <- a$Test
  ref.pvalue <- a$ref.pvalue
  margin <- a$margin
  if (margin == 0 && p1 %==% p2) stop('p1 and p2 must differ for a sample size to exist')
  if (p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1) stop('p1 and p2 must lie in [0, 1]')
  if (margin == 0 && alternative == 'greater' && p1 < p2) {
    stop("p1 must exceed p2 when alternative is 'greater'")
  }
  if (margin == 0 && alternative == 'less' && p1 > p2) {
    stop("p1 must fall below p2 when alternative is 'less'")
  }
  if (margin != 0 && alternative == 'greater' && !((p1 - p2 + margin) %>>% 0)) {
    stop("p1 - p2 must exceed -margin when alternative is 'greater'")
  }
  if (margin != 0 && alternative == 'less' && !((p2 - p1 + margin) %>>% 0)) {
    stop("p1 - p2 must fall below margin when alternative is 'less'")
  }
  if (method == 'exact' && rounding != 'group') {
    stop("rounding must be 'group' when method is 'exact'")
  }
  n <- sample_size_n(p1, p2, r, alpha, tar.power, Test, alternative, tsmethod, n.grid,
                     bb.gamma, method, rounding, ref.pvalue, margin)
  N1 <- n[['N1']]
  N2 <- n[['N2']]
  N <- N1 + N2
  rr <- get_rr(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma, ref.pvalue,
               margin)
  Power <- power_from_rr(rr, dbinom(0:N1, N1, p1), dbinom(0:N2, N2, p2))
  out <- data.frame(
    p1 = p1, p2 = p2, r = r, alpha = alpha, tar.power = tar.power,
    Test = Test, alternative = alternative, Power = Power,
    N1 = N1, N2 = N2, N = N,
    stringsAsFactors = FALSE
  )
  attr(out, 'tsmethod') <- tsmethod
  attr(out, 'n.grid') <- n.grid
  attr(out, 'bb.gamma') <- bb.gamma
  attr(out, 'method') <- method
  attr(out, 'rounding') <- rounding
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  class(out) <- c('bbssr_samplesize', 'data.frame')
  out
}
