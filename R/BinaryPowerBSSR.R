#' Power of a Blinded Sample Size Re-estimation Design for Binary Endpoints
#'
#' Calculates the power of a two-arm trial with a binary endpoint when blinded sample size
#' re-estimation (BSSR) is implemented, together with the power of the corresponding
#' fixed-sample design, the expected final sample size and its distribution. Seven tests
#' are supported (see \code{\link{BinaryRR}}), each of which can be applied with a
#' one-sided or a two-sided alternative, under either a restricted or an unrestricted design
#' rule.
#'
#' @param p Vector of true pooled proportions of responders from both groups
#' @param Delta.A Assumed treatment effect, on the scale given by \code{effect}, used to
#'   split the blinded pooled proportion into group-specific proportions
#' @param Delta.T True treatment effect, on the scale given by \code{effect}
#' @param N1 Initial sample size of group 1
#' @param N2 Initial sample size of group 2
#' @param omega Fraction of the initial sample size observed at the interim analysis. The
#'   interim size of group 2 is \code{ceiling(omega N2)} and that of group 1 is
#'   \code{ceiling(r ceiling(omega N2))}, so the interim analysis keeps the allocation
#'   ratio. Either \code{omega} or \code{n.interim} must be supplied
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance of the final analysis, for the alternative specified
#'   by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test of the final analysis. Options: \code{'Chisq'},
#'   \code{'Fisher'}, \code{'Fisher-midP'}, \code{'Z-pool'}, \code{'Boschloo'},
#'   \code{'Blackwelder'} or \code{'Farrington-Manning'}
#' @param restricted Logical. If \code{TRUE}, the re-estimated sample size is not allowed
#'   to fall below the initial sample size. Default is \code{FALSE}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default), \code{'less'} or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests, see \code{\link{BinaryRR}}. Options: \code{'minlike'} (default),
#'   \code{'central'} or \code{'blaker'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#' @param effect Scale of \code{Delta.A} and \code{Delta.T}. Options: \code{'RD'}
#'   (default) for the risk difference \code{p1 - p2}, \code{'RR'} for the risk ratio
#'   \code{p1 / p2} or \code{'OR'} for the odds ratio
#' @param ss.method How the sample size is re-estimated from the recovered proportions.
#'   Options: \code{'exact'} (default), \code{'standard'}, \code{'null.variance'} or
#'   \code{'alternative.variance'}, as in \code{\link{BinarySampleSize}}
#' @param ss.Test Test whose exact power is used for the re-estimation when
#'   \code{ss.method = 'exact'}. Default is \code{Test}
#' @param ss.alpha Level of significance used for the re-estimation. Default is
#'   \code{alpha}
#' @param rounding How the re-estimated sample size is turned into whole numbers.
#'   Options: \code{'group'} (default), \code{'friede-kieser'}, \code{'total'} or
#'   \code{'nearest'}. Only \code{'group'} is available with \code{ss.method = 'exact'}.
#'   See Details
#' @param N.min Lower bound on the final total sample size, or \code{NULL} (default) for
#'   none beyond the interim total. It can be used to keep the patients who are already
#'   enrolled but not yet evaluated at the interim analysis
#' @param N.max Upper bound on the final total sample size, or \code{NULL} (default) for
#'   none
#' @param n.interim Interim sample sizes of group 1 and group 2, as a vector of length
#'   two. An alternative to \code{omega}
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
#'   \code{\link{BinaryRR}}. Default is \code{FALSE}. It applies to the final analysis, to
#'   the fixed-sample comparator and to the exact re-estimation
#'
#' @return An object of class \code{bbssr_powerbssr}, a data frame with one row per element
#'   of \code{p} containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{p}{True pooled probability of responders from both groups}
#'   \item{power.BSSR}{Power of the BSSR design}
#'   \item{power.TRAD}{Power of the fixed-sample design}
#'   \item{E.N}{Expected total sample size of the BSSR design}
#' }
#' The interim sample sizes are stored as the attributes \code{n1.interim} and
#' \code{n2.interim}. The attribute \code{reestimation} holds, for every pooled number of
#' interim responders \code{s}, the recovered proportions and the final sample sizes, and
#' the attribute \code{N.dist} holds the distribution of the final sample size for every
#' row of the result, identified by the column \code{scenario}. \code{summary()} reports the standard deviation and quantiles of
#' the final sample size.
#'
#' @details
#' At the interim analysis the pooled number of responders is observed without unblinding.
#' The pooled proportion is combined with the assumed treatment effect \code{Delta.A} to
#' recover group-specific proportions, from which the sample size is re-estimated. For the
#' risk difference the recovered proportions are \code{hat.p + Delta.A / (1 + r)} and
#' \code{hat.p - r Delta.A / (1 + r)}. For the risk ratio they are the proportions with
#' the pooled value \code{hat.p} and the ratio \code{Delta.A}, formula (21.11) of Kieser
#' (2020), and for the odds ratio they are obtained from the quadratic equation that the
#' odds ratio and the pooled value define. The exact search of
#' \code{ss.method = 'exact'} uses the recovered proportions truncated to the unit
#' interval. The normal approximation uses them before truncation, with each Bernoulli
#' variance truncated at zero, so that the assumed effect is kept as in formula (2) of
#' Friede and Kieser (2004). The power is then averaged over the distribution of the
#' interim outcome.
#'
#' Under \code{rounding = 'group'}, both the interim and the final sample size of group 1
#' are obtained from the size of group 2 by a single application of
#' \code{ceiling(r ...)}, so the allocation ratio is preserved as closely as whole numbers
#' allow and is exact whenever \code{r} is a whole number. Under
#' \code{rounding = 'friede-kieser'}, which requires a normal approximation, that is an
#' \code{ss.method} other than \code{'exact'}, the second stage receives the rounded-up
#' excess of the unrounded re-estimated total over the interim total, split into
#' \code{ceiling(n / (1 + r))} and \code{ceiling(r n / (1 + r))} patients as in Friede and
#' Kieser (2004). Under \code{rounding = 'total'} the unrounded total is rounded up and
#' group 2 receives \code{floor(N / (1 + r))} patients. Under \code{rounding = 'nearest'}
#' the unrounded total is split in the ratio \code{r} to 1 and each group is rounded to
#' the nearest whole number, which reproduces the computations of Friede et al. (2007).
#' In each case the final total is kept between the interim total (or \code{N.min}, or
#' the initial total under the restricted rule) and \code{N.max}, except that
#' \code{'friede-kieser'} and \code{'nearest'} round the two groups separately and can
#' exceed \code{N.max} by one patient. The argument
#' \code{N1} enters only through the
#' fixed-sample comparator and the restricted rule, and under \code{rounding = 'group'} a
#' warning is issued when it is not \code{ceiling(r N2)}.
#'
#' Recovered proportions that coincide admit no sample size. The risk ratio produces them
#' when no interim patient responds, and the odds ratio when no interim patient or every
#' interim patient responds. The initial sample size is kept for such interim outcomes.
#'
#' Setting \code{Delta.T} to the null value (0 for the risk difference, 1 for a ratio)
#' makes the two groups identical, so \code{power.BSSR} and \code{power.TRAD} become
#' rejection probabilities under the null hypothesis. \code{\link{BinaryTypeIErrorBSSR}}
#' evaluates them over the unit interval and locates the maximum.
#'
#' With a non-inferiority \code{margin}, the design of Friede et al. (2007) is obtained.
#' The assumed effect \code{Delta.A} is usually 0, which recovers the pooled proportion
#' for both groups, the sample size is re-estimated by the formula of Farrington and
#' Manning (\code{ss.method = 'standard'}) or of Blackwelder
#' (\code{ss.method = 'alternative.variance'}), and \code{Delta.T = -margin} (or
#' \code{margin} for \code{alternative = 'less'}) gives the rejection probability on the
#' boundary of the null hypothesis.
#'
#' @references
#' Friede T, Kieser M (2004). Sample size recalculation for binary data in internal pilot
#' study designs. \emph{Pharmaceutical Statistics}, 3(4), 269-279.
#'
#' Kieser M (2020). \emph{Methods and Applications of Sample Size Calculation and
#' Recalculation in Clinical Trials}. Springer, Cham.
#'
#' Friede T, Mitchell C, Mueller-Velten G (2007). Blinded sample size reestimation in
#' non-inferiority trials with binary endpoints. \emph{Biometrical Journal}, 49(6),
#' 903-916.
#'
#' @examples
#' # Small BSSR calculation with the chi-squared test
#' BinaryPowerBSSR(
#'   p = 0.45,
#'   Delta.A = 0.3, Delta.T = 0.3,
#'   N1 = 5, N2 = 5, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#'
#' \donttest{
#' res <- BinaryPowerBSSR(
#'   p = seq(0.19, 0.37, by = 0.03),
#'   Delta.A = 0.36, Delta.T = 0.36,
#'   N1 = 24, N2 = 24, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'
#' )
#' print(res)
#' summary(res)
#' plot(res)
#'
#' # Re-estimation by the normal approximation with an upper bound of twice the
#' # initial sample size, and a final analysis with the Boschloo test
#' BinaryPowerBSSR(
#'   p = seq(0.2, 0.4, by = 0.05),
#'   Delta.A = 0.3, Delta.T = 0.3,
#'   N1 = 30, N2 = 30, n.interim = c(15, 15), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Boschloo',
#'   ss.method = 'standard', N.max = 120
#' )
#'
#' # Non-inferiority with a margin of 0.1, as in Friede et al. (2007)
#' BinaryPowerBSSR(
#'   p = c(0.5, 0.7), Delta.A = 0, Delta.T = 0,
#'   N1 = 329, N2 = 329, n.interim = c(66, 66), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Farrington-Manning',
#'   ss.method = 'standard', rounding = 'nearest', margin = 0.1
#' )
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryTypeIErrorBSSR}}, \code{\link{BinaryAlphaAdjBSSR}},
#'   \code{\link{BinaryBSSR}}
#' @export
#' @import fpCompare
#' @importFrom stats dbinom
BinaryPowerBSSR <- function(p, Delta.A, Delta.T, N1, N2, omega = NULL, r,
                            alpha, tar.power, Test, restricted = FALSE,
                            alternative = c('greater', 'less', 'two.sided'),
                            tsmethod = c('minlike', 'central', 'blaker'),
                            n.grid = 100, bb.gamma = 0,
                            effect = c('RD', 'RR', 'OR'),
                            ss.method = c('exact', 'standard', 'null.variance',
                                          'alternative.variance'),
                            ss.Test = Test, ss.alpha = alpha,
                            rounding = c('group', 'friede-kieser', 'total', 'nearest'),
                            N.min = NULL, N.max = NULL, n.interim = NULL, margin = 0,
                            ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  if (rounding == 'group' && N1 != ceiling(r * N2)) {
    warning('N1 differs from ceiling(r * N2), the fixed-sample comparator is not ',
            'allocated in the ratio r to 1')
  }
  # Final sample size for every pooled number of interim responders
  map <- bssr_map(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                  restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                  ss.method, ss.Test, ss.alpha, rounding, N.min, N.max, ref.pvalue,
                  margin)
  setup <- bssr_setup(map)
  N11 <- setup$n11
  N12 <- setup$n12
  # True proportion for each treatment group
  sp <- split_pooled(p, Delta.T, r, effect)
  p1 <- sp$p1
  p2 <- sp$p2
  # Omit scenarios with p_{j} outside the unit interval
  omit.ID <- which(p1 %>>% 1 | p1 %<<% 0 | p2 %>>% 1 | p2 %<<% 0)
  if (length(omit.ID) %!=% 0) {
    p1 <- p1[-omit.ID]
    p2 <- p2[-omit.ID]
    p <- p[-omit.ID]
  }
  if (length(p) == 0) stop('no scenario has both p1 and p2 inside the unit interval')
  # Probabilities that lie outside the unit interval only by rounding error are kept, and
  # are moved onto it so that the binomial probabilities are defined
  p1 <- pmin(1, pmax(0, p1))
  p2 <- pmin(1, pmax(0, p2))
  # Probability mass functions of the interim outcome
  dbinom1 <- outer(X = 0:N11, Y = p1, function(X, Y) dbinom(X, N11, Y))
  dbinom2 <- outer(X = 0:N12, Y = p2, function(X, Y) dbinom(X, N12, Y))
  # Rejection regions of the distinct final sample sizes, and the power of the design.
  # The conditional power of the second stage is summed in compiled code over the runs of
  # rejected outcomes in each column of the rejection region
  rr.list <- lapply(seq_along(setup$N1), function(k) {
    get_rr(setup$N1[k], setup$N2[k], alpha, Test, alternative, tsmethod, n.grid, bb.gamma,
           ref.pvalue, margin)
  })
  power.BSSR <- bssr_reject(setup, rr.list, p1, p2)
  # Final total sample size of every interim outcome, in the order of the interim cells
  s.cell <- setup$x11 + setup$x12
  hat.N <- (map$N1 + map$N2)[s.cell + 1L]
  E.N <- vapply(seq_along(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * hat.N)
  }, numeric(1))
  # Distribution of the final sample size
  N.dist <- do.call(rbind, lapply(seq_along(p), function(k) {
    prob.s <- rowsum(c(dbinom1[, k] %o% dbinom2[, k]), s.cell, reorder = TRUE)[, 1]
    key <- paste(map$N1, map$N2, sep = '_')
    prob.N <- rowsum(prob.s, key, reorder = FALSE)[, 1]
    first <- match(names(prob.N), key)
    d <- data.frame(scenario = k, p = p[k], N1 = map$N1[first], N2 = map$N2[first],
                    N = map$N1[first] + map$N2[first], prob = unname(prob.N))
    d[order(d$N, d$N1), , drop = FALSE]
  }))
  rownames(N.dist) <- NULL
  # Power of the fixed-sample design
  power.TRAD <- BinaryPower(p1, p2, N1, N2, alpha, Test, alternative, tsmethod,
                            n.grid, bb.gamma, margin = margin,
                            ref.pvalue = ref.pvalue)$Power
  out <- data.frame(p1, p2, p, power.BSSR, power.TRAD, E.N)
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N1') <- N1
  attr(out, 'N2') <- N2
  attr(out, 'omega') <- omega
  attr(out, 'n1.interim') <- N11
  attr(out, 'n2.interim') <- N12
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'Delta.T') <- Delta.T
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'rounding') <- rounding
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  attr(out, 'reestimation') <- map
  attr(out, 'N.dist') <- N.dist
  class(out) <- c('bbssr_powerbssr', 'data.frame')
  out
}
