#' Blinded Sample Size Re-Estimation from Paired Interim Data
#'
#' Re-estimates the number of pairs of a matched-pairs trial from the number of discordant
#' pairs observed at an interim analysis, without the direction in which they are
#' discordant.
#'
#' @param n1 Number of pairs at the interim analysis
#' @param D1 Number of discordant pairs among them
#' @param N Planned number of pairs. Required when \code{restricted = TRUE}, otherwise
#'   optional
#' @inheritParams PairedPowerBSSR
#'
#' @return An object of class \code{bbssr_pbssr}, a data frame with one row containing
#'   the interim data \code{n1} and \code{D1}, the blinded estimate \code{psi.hat} of the
#'   probability of a discordant pair, the value \code{psi.plan} used in the
#'   re-estimation, the planning probabilities \code{p10.plan} and \code{p01.plan}, the
#'   re-estimated number of pairs \code{N.re} before the bounds are applied (at most
#'   \code{N.max} with \code{ss.method = 'exact'}, whose search does not go beyond it, and
#'   within the bounds with \code{ss.method = 'conditional.power'}), the number
#'   \code{n.stage2} of pairs still to be enrolled, the final number of pairs
#'   \code{N.final}, the power \code{Power} of the final analysis with \code{N.final}
#'   pairs under the planning probabilities, and its conditional power \code{CP} given
#'   \code{D1}, computed with \code{ss.Test} at the level \code{ss.alpha}. The attributes
#'   hold the design settings.
#'
#' @details
#' The re-estimation is that of \code{\link{PairedPowerBSSR}}: the blinded estimate
#' \code{D1 / n1} of the probability of a discordant pair, raised to \code{|Delta.A|} when
#' it is smaller, and the assumed difference \code{Delta.A} give the planning
#' probabilities, from which the number of pairs follows by the rule \code{ss.method}.
#' The conditional power is
#' \code{sum over D2 of b(D2; N.final - n1, psi.plan) g(D1 + D2)}, the probability of
#' rejection given \code{D1} under the planning probabilities, formula (4) of Schepers et
#' al. (2026) for the exact conditional McNemar test.
#'
#' @references
#' Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026). Blinded sample size
#' review for McNemar's test based on primary and surrogate endpoints. arXiv:2608.17784.
#'
#' @examples
#' # 25 interim pairs of which 5 are discordant, assumed difference 0.15
#' PairedBSSR(n1 = 25, D1 = 5, Delta.A = 0.15, alpha = 0.05, tar.power = 0.8,
#'            Test = 'McNemar-exact', alternative = 'two.sided',
#'            ss.method = 'conditional.power', N.max = 150)
#'
#' @export
PairedBSSR <- function(n1, D1, Delta.A, alpha, tar.power, Test, restricted = FALSE,
                       N = NULL,
                       alternative = c('greater', 'less', 'two.sided', 'equivalence'),
                       n.grid = 100, bb.gamma = 0,
                       ss.method = c('exact', 'standard', 'conditional.power'),
                       ss.Test = Test, ss.alpha = alpha,
                       search = c('crossing', 'smallest', 'stable'),
                       search.limit = c(2, 50), N.min = NULL, N.max = NULL, margin = 0,
                       ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  margin <- norm_margin(alternative, margin, 'RD')
  if (length(n1) != 1 || length(D1) != 1) stop('n1 and D1 must each be a single value')
  if (is.na(n1) || n1 != round(n1) || n1 < 1) stop('n1 must be a positive integer')
  if (is.na(D1) || D1 != round(D1) || D1 < 0 || D1 > n1) {
    stop('D1 must be an integer between 0 and n1')
  }
  if (!is.null(N) && (length(N) != 1 || is.na(N) || N < 1 || N != round(N))) {
    stop('N must be a single positive integer')
  }
  # The test of the final analysis, with the Berger-Boos parameter only for the
  # unconditional test
  Test <- check_paired_args(1, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                            margin)$Test
  bb.test <- if (Test == 'McNemar-uncond') bb.gamma else 0
  re <- paired_reestimate(as.integer(n1), as.integer(D1), Delta.A, N, tar.power,
                          restricted, alternative, n.grid, bb.gamma, ss.method, ss.Test,
                          ss.alpha, N.min, N.max, ref.pvalue, margin, search, search.limit)
  N.final <- re$N
  Power <- PairedPower(re$p10.plan, re$p01.plan, N.final, alpha, Test, alternative,
                       n.grid, bb.test, margin, ref.pvalue)$Power
  ss.Test <- check_paired_args(1, ss.alpha, ss.Test, n.grid, 0, ref.pvalue, alternative,
                               margin)$Test
  rr <- get_paired_rr(N.final, ss.alpha, ss.Test, alternative, n.grid,
                      if (ss.Test == 'McNemar-uncond') bb.gamma else 0, ref.pvalue, margin)
  pi.dir <- if (re$psi.plan > 0) re$p10.plan / re$psi.plan else 0.5
  CP <- paired_cp(rr, as.integer(n1), as.integer(D1), re$psi.plan, pi.dir)
  out <- data.frame(n1 = as.integer(n1), D1 = as.integer(D1), psi.hat = re$psi.hat,
                    psi.plan = re$psi.plan, p10.plan = re$p10.plan,
                    p01.plan = re$p01.plan, N.re = re$N.re,
                    n.stage2 = as.integer(N.final - n1), N.final = N.final,
                    Power = Power, CP = CP)
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N') <- N
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'search') <- search
  attr(out, 'search.limit') <- attr(re, 'N.limit')
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  class(out) <- c('bbssr_pbssr', 'data.frame')
  out
}
