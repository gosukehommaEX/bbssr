#' Sample Size Re-estimation from Observed Blinded Interim Data
#'
#' Re-estimates the sample size of an ongoing two-arm trial with a binary endpoint from the
#' blinded data available at an interim analysis, and reports how many patients still have
#' to be enrolled in each group during the second stage. Only the total number of patients
#' and the total number of responders are required, so the treatment allocation remains
#' concealed.
#'
#' @param n1 Number of patients of group 1 observed at the interim analysis
#' @param n2 Number of patients of group 2 observed at the interim analysis
#' @param S Total number of responders observed at the interim analysis, pooled over both
#'   groups
#' @param Delta.A Assumed treatment effect, on the scale given by \code{effect}, used to
#'   split the blinded pooled proportion into group-specific proportions
#' @param r Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2 = r:1,
#'   r > 0)
#' @param alpha Level of significance of the final analysis, for the alternative specified
#'   by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test of the final analysis. Options: \code{'Chisq'},
#'   \code{'Fisher'}, \code{'Fisher-midP'}, \code{'Z-pool'}, \code{'Boschloo'},
#'   \code{'Blackwelder'} or \code{'Farrington-Manning'}
#' @param restricted Logical. If \code{TRUE}, the re-estimated sample size is not allowed
#'   to fall below the planned sample size given by \code{N1} and \code{N2}. Default is
#'   \code{FALSE}
#' @param N1 Planned sample size of group 1. Required when \code{restricted} is \code{TRUE},
#'   and when the recovered proportions coincide
#' @param N2 Planned sample size of group 2. Required when \code{restricted} is \code{TRUE},
#'   and when the recovered proportions coincide
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default), \code{'less'} or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests, see \code{\link{BinaryRR}}. Options: \code{'minlike'} (default),
#'   \code{'central'} or \code{'blaker'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#' @param effect Scale of \code{Delta.A}. Options: \code{'RD'} (default), \code{'RR'} or
#'   \code{'OR'}, as in \code{\link{BinaryPowerBSSR}}
#' @param ss.method How the sample size is re-estimated. Options: \code{'exact'}
#'   (default), \code{'standard'}, \code{'null.variance'} or
#'   \code{'alternative.variance'}, as in \code{\link{BinarySampleSize}}
#' @param ss.Test Test whose exact power is used for the re-estimation when
#'   \code{ss.method = 'exact'}. Default is \code{Test}
#' @param ss.alpha Level of significance used for the re-estimation. Default is
#'   \code{alpha}
#' @param rounding How the re-estimated sample size is turned into whole numbers. Options:
#'   \code{'group'} (default), \code{'friede-kieser'}, \code{'total'} or \code{'nearest'},
#'   as in \code{\link{BinaryPowerBSSR}}
#' @param N.min Lower bound on the final total sample size, or \code{NULL} (default). It
#'   can be used to keep the patients who are already enrolled but not yet evaluated
#' @param N.max Upper bound on the final total sample size, or \code{NULL} (default)
#' @param margin Non-inferiority margin on the scale of the risk difference. The
#'   default of 0 gives a test of superiority. A value other than 0 tests the null
#'   hypothesis \code{p1 - p2 <= -margin} against \code{p1 - p2 > -margin} when
#'   \code{alternative} is \code{'greater'}, and \code{p1 - p2 >= margin} against
#'   \code{p1 - p2 < margin} when it is \code{'less'}. It requires
#'   \code{Test = 'Blackwelder'} or \code{'Farrington-Manning'}, see
#'   \code{\link{BinaryRR}}. A negative value tests for superiority by more than its
#'   absolute value. With a value other than 0 the assumed
#'   and the true effects are risk differences (\code{effect = 'RD'})
#' @param ref.pvalue Logical. If \code{TRUE}, the maximization over the nuisance parameter
#'   of the unconditional tests is refined between the grid points, see
#'   \code{\link{BinaryRR}}. Default is \code{FALSE}
#'
#' @return An object of class \code{bbssr_bssr}, a data frame with one row containing:
#' \describe{
#'   \item{n1}{Interim sample size of group 1}
#'   \item{n2}{Interim sample size of group 2}
#'   \item{n}{Total interim sample size}
#'   \item{S}{Total number of interim responders}
#'   \item{hat.p}{Blinded estimate of the pooled response probability}
#'   \item{hat.p1}{Recovered response probability of group 1}
#'   \item{hat.p2}{Recovered response probability of group 2}
#'   \item{N1.re}{Re-estimated total sample size of group 1}
#'   \item{N2.re}{Re-estimated total sample size of group 2}
#'   \item{N.re}{Re-estimated total sample size}
#'   \item{n1.stage2}{Number of additional patients to enrol in group 1}
#'   \item{n2.stage2}{Number of additional patients to enrol in group 2}
#'   \item{n.stage2}{Total number of additional patients to enrol}
#'   \item{N1.final}{Final sample size of group 1}
#'   \item{N2.final}{Final sample size of group 2}
#'   \item{N.final}{Final total sample size}
#'   \item{Power}{Exact power at the final sample size under the recovered proportions}
#' }
#'
#' @details
#' The blinded estimate of the pooled response probability is \code{hat.p = S / (n1 + n2)}.
#' For the risk difference, group-specific proportions are recovered as
#' \code{hat.p1 = hat.p + Delta.A / (1 + r)} and
#' \code{hat.p2 = hat.p - r Delta.A / (1 + r)}, truncated to the unit interval. For the
#' risk ratio and the odds ratio they are the proportions with the pooled value
#' \code{hat.p} and the ratio \code{Delta.A}. The sample size is then re-estimated as in
#' \code{\link{BinarySampleSize}}.
#'
#' Under the unrestricted rule and \code{rounding = 'group'}, the final size of group 2 is
#' the larger of the re-estimated size and what has already been observed. Under the
#' restricted rule it is raised to the planned size first, so the trial can only grow, and
#' \code{N.min} and \code{N.max} add a further lower and upper bound on the total. The
#' final size of group 1 is then \code{ceiling(r N2.final)}, so an imbalance already
#' present at the interim is corrected by the remaining enrolment instead of being carried
#' forward. Neither second-stage size is ever negative. The other rounding rules are
#' described in \code{\link{BinaryPowerBSSR}}.
#'
#' While \code{\link{BinaryPowerBSSR}} evaluates the operating characteristics of a BSSR
#' design at the planning stage, this function is applied once, to the data of a trial that
#' is under way.
#'
#' @examples
#' # Interim data: 20 patients per group, 11 responders in total
#' BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
#'            alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
#'
#' \donttest{
#' # Restricted rule with a planned sample size of 40 per group
#' BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
#'            alpha = 0.025, tar.power = 0.8, Test = 'Boschloo',
#'            restricted = TRUE, N1 = 40, N2 = 40)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryPowerBSSR}}, \code{\link{BinarySampleSize}}
#' @export
BinaryBSSR <- function(n1, n2, S, Delta.A, r, alpha, tar.power, Test,
                       restricted = FALSE, N1 = NULL, N2 = NULL,
                       alternative = c('greater', 'less', 'two.sided'),
                       tsmethod = c('minlike', 'central', 'blaker'),
                       n.grid = 100, bb.gamma = 0,
                       effect = c('RD', 'RR', 'OR'),
                       ss.method = c('exact', 'standard', 'null.variance',
                                     'alternative.variance'),
                       ss.Test = Test, ss.alpha = alpha,
                       rounding = c('group', 'friede-kieser', 'total', 'nearest'),
                       N.min = NULL, N.max = NULL, margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  if (length(n1) != 1 || length(n2) != 1 || length(S) != 1) {
    stop('n1, n2 and S must each be a single value')
  }
  if (n1 != round(n1) || n2 != round(n2) || n1 < 1 || n2 < 1) {
    stop('n1 and n2 must be positive integers')
  }
  if (S != round(S) || S < 0 || S > n1 + n2) {
    stop('S must be an integer between 0 and n1 + n2')
  }
  check_delta(Delta.A, effect, alternative, 'Delta.A', margin)
  if (restricted && (is.null(N1) || is.null(N2))) {
    stop('N1 and N2 must be supplied when restricted is TRUE')
  }
  if (ss.method == 'exact' && rounding != 'group') {
    stop("rounding must be 'group' when ss.method is 'exact'")
  }
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  S <- as.integer(S)
  # Blinded estimate of the pooled response probability
  hat.p <- S / (n1 + n2)
  sp <- split_pooled(hat.p, Delta.A, r, effect)
  hat.p1 <- pmin(1, pmax(0, sp$p1))
  hat.p2 <- pmin(1, pmax(0, sp$p2))
  # Sample size re-estimation
  # The exact search needs probabilities inside the unit interval. The normal
  # approximation uses the untruncated probabilities, which keep the assumed effect
  ss.p <- if (ss.method == 'exact') list(p1 = hat.p1, p2 = hat.p2) else sp
  re <- reestimate(ss.p$p1, ss.p$p2, r, ss.alpha, tar.power, ss.Test, alternative,
                   tsmethod, n.grid, bb.gamma, ss.method, rounding, N1, N2, ref.pvalue,
                   margin)
  N1.re <- re$N1.re
  N2.re <- re$N2.re
  # Final sample sizes. Under the group rounding the final size of group 2 is fixed
  # first, and group 1 is brought to ceiling(r N2.final), so any imbalance already present
  # at the interim is corrected by the remaining enrolment rather than carried forward
  fin <- final_sizes(re, n1, n2, r, rounding, restricted, N1, N2, N.min, N.max)
  N1.final <- fin$N1
  N2.final <- fin$N2
  n1.stage2 <- N1.final - n1
  n2.stage2 <- N2.final - n2
  Power <- BinaryPower(hat.p1, hat.p2, N1.final, N2.final, alpha, Test,
                       alternative, tsmethod, n.grid, bb.gamma, margin = margin,
                       ref.pvalue = ref.pvalue)$Power
  out <- data.frame(
    n1 = n1, n2 = n2, n = n1 + n2, S = S,
    hat.p = hat.p, hat.p1 = hat.p1, hat.p2 = hat.p2,
    N1.re = N1.re, N2.re = N2.re, N.re = N1.re + N2.re,
    n1.stage2 = n1.stage2, n2.stage2 = n2.stage2, n.stage2 = n1.stage2 + n2.stage2,
    N1.final = N1.final, N2.final = N2.final, N.final = N1.final + N2.final,
    Power = Power,
    stringsAsFactors = FALSE
  )
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'r') <- r
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'rounding') <- rounding
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  class(out) <- c('bbssr_bssr', 'data.frame')
  out
}
