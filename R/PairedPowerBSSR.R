#' Power of Blinded Sample Size Re-Estimation for Paired Binary Data
#'
#' Calculates the exact power and the distribution of the final number of pairs of a
#' matched-pairs design with blinded sample size re-estimation, together with the power of
#' the fixed design with the planned number of pairs. The number of pairs is re-estimated
#' at an interim analysis from the number of discordant pairs, without the direction in
#' which they are discordant, and the final analysis uses all pairs.
#'
#' @param psi True probability \code{p10 + p01} of a discordant pair. A vector is allowed
#' @param Delta.T True difference \code{p10 - p01}, a single value
#' @param Delta.A Assumed difference \code{p10 - p01} used in the re-estimation
#' @param N Planned number of pairs, which is the number of pairs of the fixed design
#' @param omega Fraction of \code{N} at the interim analysis, rounded up. Supply either
#'   \code{omega} or \code{n.interim}
#' @param tar.power Target power of the re-estimation
#' @param restricted Logical. If \code{TRUE}, the final number of pairs is at least
#'   \code{N}. Default is \code{FALSE}
#' @param ss.method Re-estimation rule. \code{'exact'} (default) and \code{'standard'}
#'   re-estimate the number of pairs of \code{\link{PairedSampleSize}} with the method of
#'   that name. \code{'conditional.power'} takes the smallest number of pairs at which the
#'   conditional power given the interim number of discordant pairs attains
#'   \code{tar.power}, the rule of Schepers et al. (2026); it requires \code{N.max}. See
#'   Details
#' @param ss.Test Test used in the re-estimation. Default is \code{Test}
#' @param ss.alpha Level used in the re-estimation. Default is \code{alpha}
#' @param search Search of the exact re-estimation, see \code{\link{PairedSampleSize}};
#'   it does not go beyond \code{N.max}
#' @param search.limit Limit of the search \code{'stable'}, see
#'   \code{\link{PairedSampleSize}}
#' @param N.min Smallest final number of pairs, or \code{NULL}
#' @param N.max Largest final number of pairs, or \code{NULL}
#' @param n.interim Number of pairs at the interim analysis. Supply either \code{omega} or
#'   \code{n.interim}
#' @inheritParams PairedRR
#'
#' @return An object of class \code{bbssr_ppowerbssr}, a data frame with one row for each
#'   element of \code{psi} that gives probabilities \code{p10} and \code{p01} in
#'   \code{[0, 1]}, containing:
#' \describe{
#'   \item{psi, p10, p01}{True probabilities of a discordant pair and of its two kinds}
#'   \item{power.BSSR}{Power of the design with re-estimation}
#'   \item{power.TRAD}{Power of the fixed design with \code{N} pairs}
#'   \item{E.N}{Expected final number of pairs}
#' }
#' The attribute \code{reestimation} holds the final number of pairs for every interim
#' number of discordant pairs, \code{N.dist} the distribution of the final number of pairs
#' in each scenario, and the other attributes the design settings.
#'
#' @details
#' Among the \code{n1} pairs of the interim analysis, the number \code{D1} of discordant
#' pairs is binomial with \code{n1} trials and probability \code{psi}, and the blinded
#' estimate of \code{psi} is \code{D1 / n1}. The assumed difference \code{Delta.A} is
#' kept, and the estimate is raised to \code{|Delta.A|} when it is smaller, since
#' \code{psi >= |p10 - p01|}. The planning probabilities are then
#' \code{p10 = (psi.plan + Delta.A) / 2} and \code{p01 = (psi.plan - Delta.A) / 2}, as in
#' Xu and Yu (2013) and Schepers et al. (2026).
#'
#' The final number of pairs depends on the interim data only through \code{D1}. Given the
#' numbers of discordant pairs of the two stages, the number of pairs in which only
#' treatment 1 responds is binomial with the total number of discordant pairs and
#' probability \code{p10 / psi}, whatever the re-estimation did. The power is therefore
#' \code{sum over D1 and D2 of b(D1; n1, psi) b(D2; N(D1) - n1, psi) g(D1 + D2)}, where
#' \code{g(D)} is the probability that the final analysis with \code{N(D1)} pairs rejects
#' given \code{D} discordant pairs. It is computed exactly by enumeration.
#'
#' With \code{ss.method = 'conditional.power'} the final number of pairs is the smallest
#' number \code{N} between the bounds at which
#' \code{sum over D2 of b(D2; N - n1, psi.plan) g_N(D1 + D2)}, computed with the planning
#' probabilities and the test \code{ss.Test} at the level \code{ss.alpha}, attains
#' \code{tar.power}. With the exact conditional McNemar test this is the conditional power
#' (4) of Schepers et al. (2026), which they evaluate with \code{p* = 0.5 + Delta / (2 psi)}
#' at the blinded estimate of \code{psi}. When no number up to \code{N.max} attains the
#' target, \code{N.max} is used, a choice the article leaves open.
#'
#' The exact conditional McNemar test keeps its level whatever blinded rule is used,
#' since the number of pairs in which only treatment 1 responds stays binomial given the
#' total number of discordant pairs (Schepers et al., 2026). The other tests can exceed
#' their nominal level after re-estimation, see \code{\link{PairedTypeIErrorBSSR}}.
#'
#' @references
#' Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026). Blinded sample size
#' review for McNemar's test based on primary and surrogate endpoints. arXiv:2608.17784.
#'
#' Xu J, Yu M (2013). Sample size determination and re-estimation for matched pair designs
#' with multiple binary endpoints. \emph{Biometrical Journal}, 55(3), 430-443.
#'
#' @examples
#' res <- PairedPowerBSSR(psi = c(0.15, 0.2, 0.3), Delta.T = 0.15, Delta.A = 0.15, N = 72,
#'                        n.interim = 30, alpha = 0.05, tar.power = 0.8,
#'                        Test = 'McNemar-exact', alternative = 'two.sided',
#'                        ss.method = 'standard')
#' print(res)
#'
#' \donttest{
#' # Conditional power rule of Schepers et al. (2026) with at most 150 pairs
#' res <- PairedPowerBSSR(psi = c(0.15, 0.2, 0.3), Delta.T = 0.15, Delta.A = 0.15, N = 72,
#'                        n.interim = 30, alpha = 0.05, tar.power = 0.8,
#'                        Test = 'McNemar-exact', alternative = 'two.sided',
#'                        ss.method = 'conditional.power', N.max = 150)
#' print(res)
#' }
#'
#' @export
#' @importFrom stats dbinom
PairedPowerBSSR <- function(psi, Delta.T, Delta.A, N, omega = NULL, alpha, tar.power,
                            Test, restricted = FALSE,
                            alternative = c('greater', 'less', 'two.sided', 'equivalence'),
                            n.grid = 100, bb.gamma = 0,
                            ss.method = c('exact', 'standard', 'conditional.power'),
                            ss.Test = Test, ss.alpha = alpha,
                            search = c('crossing', 'smallest', 'stable'),
                            search.limit = c(2, 50), N.min = NULL, N.max = NULL,
                            n.interim = NULL, margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  margin <- norm_margin(alternative, margin, 'RD')
  if (length(N) != 1 || is.na(N) || N < 1 || N != round(N)) {
    stop('N must be a single positive integer')
  }
  if (length(Delta.T) != 1 || is.na(Delta.T) || abs(Delta.T) >= 1) {
    stop('Delta.T must be a single value in (-1, 1)')
  }
  if (!is.numeric(psi) || length(psi) == 0 || anyNA(psi)) {
    stop('psi must be a numeric vector without missing values')
  }
  n1 <- paired_n1(N, omega, n.interim)
  # The test of the final analysis, with the Berger-Boos parameter only for the
  # unconditional test
  Test <- check_paired_args(N, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                            margin)$Test
  bb.test <- if (Test == 'McNemar-uncond') bb.gamma else 0
  # Final number of pairs for every interim number of discordant pairs
  map <- paired_reestimate(n1, 0:n1, Delta.A, N, tar.power, restricted, alternative,
                           n.grid, bb.gamma, ss.method, ss.Test, ss.alpha, N.min, N.max,
                           ref.pvalue, margin, search, search.limit)
  Ns <- sort(unique(map$N))
  rr.id <- as.integer(match(map$N, Ns) - 1L)
  rr.list <- lapply(Ns, function(n) {
    get_paired_rr(n, alpha, Test, alternative, n.grid, bb.test, ref.pvalue, margin)
  })
  # True probabilities of the two kinds of discordant pairs. Scenarios with psi below
  # |Delta.T| or above 1 are omitted
  p10 <- (psi + Delta.T) / 2
  p01 <- (psi - Delta.T) / 2
  keep <- !(p10 %<<% 0 | p01 %<<% 0 | psi %>>% 1)
  psi <- psi[keep]
  if (length(psi) == 0) {
    stop('no value of psi gives probabilities p10 and p01 in [0, 1]; psi must lie ',
         'between |Delta.T| and 1')
  }
  p10 <- pmax(0, p10[keep])
  p01 <- pmax(0, p01[keep])
  power.BSSR <- paired_bssr_power(rr.list, rr.id, n1, p10, p01)
  # Distribution of the interim number of discordant pairs in each scenario
  dD1 <- outer(0:n1, psi, function(x, t) dbinom(x, n1, pmin(1, t)))
  E.N <- colSums(dD1 * map$N)
  N.dist <- do.call(rbind, lapply(seq_along(psi), function(k) {
    prob <- rowsum(dD1[, k], map$N, reorder = TRUE)[, 1]
    data.frame(scenario = k, psi = psi[k], N = as.integer(names(prob)),
               prob = unname(prob))
  }))
  rownames(N.dist) <- NULL
  power.TRAD <- PairedPower(p10, p01, N, alpha, Test, alternative, n.grid, bb.test,
                            margin, ref.pvalue)$Power
  out <- data.frame(psi = psi, p10 = p10, p01 = p01, power.BSSR = power.BSSR,
                    power.TRAD = power.TRAD, E.N = E.N)
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N') <- as.integer(N)
  attr(out, 'omega') <- omega
  attr(out, 'n.interim') <- n1
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'Delta.T') <- Delta.T
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'search') <- search
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  attr(out, 'reestimation') <- map
  attr(out, 'N.dist') <- N.dist
  class(out) <- c('bbssr_ppowerbssr', 'data.frame')
  out
}
