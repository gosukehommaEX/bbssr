#' Re-Estimated Number of Pairs from Blinded Interim Data
#'
#' Internal helper returning the final number of pairs of a design with blinded sample
#' size re-estimation for paired binary data, for each number \code{D1} of discordant
#' pairs among the \code{n1} pairs of the interim analysis. The blinded estimate of the
#' probability of a discordant pair is \code{psi.hat = D1 / n1}. The assumed difference
#' \code{Delta.A = p10 - p01} is kept, and since \code{psi >= |p10 - p01|} the estimate
#' is raised to \code{psi.plan = max(psi.hat, |Delta.A|)}, which gives the planning
#' probabilities \code{p10 = (psi.plan + Delta.A) / 2} and
#' \code{p01 = (psi.plan - Delta.A) / 2}.
#'
#' With \code{ss.method = 'exact'} or \code{'standard'} the number of pairs is that of
#' \code{\link{PairedSampleSize}} for the planning probabilities, with the test
#' \code{ss.Test} at the level \code{ss.alpha}. With \code{ss.method = 'conditional.power'}
#' it is the smallest number \code{N} of pairs in the allowed range at which the
#' conditional power of \code{ss.Test} at the level \code{ss.alpha}, given \code{D1} and
#' the planning probabilities, attains \code{tar.power}, the rule of Schepers et al.
#' (2026). The conditional power is
#' \code{sum over D2 of b(D2; N - n1, psi.plan) g_N(D1 + D2)}, where \code{g_N(D)} is the
#' probability of rejection given \code{D} discordant pairs among the \code{N} pairs, with
#' the number of pairs in which only treatment 1 responds binomial with \code{D} trials and
#' probability \code{p10 / psi.plan}. It depends on the interim data only through
#' \code{D1}, and for the exact conditional McNemar test it is formula (4) of Schepers et
#' al. (2026). When no number up to the upper bound attains the target, the upper bound is
#' used.
#'
#' The final number of pairs is at least \code{n1}, at least \code{N} when
#' \code{restricted = TRUE} and at least \code{N.min}, and at most \code{N.max}.
#'
#' @param n1 Number of pairs at the interim analysis
#' @param D1 Numbers of discordant pairs at the interim analysis, a vector
#' @inheritParams PairedPowerBSSR
#'
#' @return A data frame with one row for each element of \code{D1}, containing \code{D1},
#'   \code{psi.hat}, \code{psi.plan}, \code{p10.plan}, \code{p01.plan}, the unrounded
#'   normal approximation \code{n.raw} (\code{NA} unless \code{ss.method = 'standard'}),
#'   the re-estimated number of pairs \code{N.re} before the bounds are applied (with
#'   \code{ss.method = 'conditional.power'} the search is confined to the bounds, so
#'   \code{N.re} already lies within them), and the final number of pairs \code{N}. The
#'   attribute \code{N.limit} holds the limits of the search \code{'stable'}
#'
#' @references
#' Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026). Blinded sample size
#' review for McNemar's test based on primary and surrogate endpoints. arXiv:2608.17784.
#'
#' @keywords internal
#' @noRd
paired_reestimate <- function(n1, D1, Delta.A, N, tar.power, restricted, alternative,
                              n.grid, bb.gamma, ss.method, ss.Test, ss.alpha, N.min,
                              N.max, ref.pvalue, margin, search, search.limit) {
  if (length(tar.power) != 1 || is.na(tar.power) || tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  if (!is.null(N.min) && (length(N.min) != 1 || is.na(N.min))) {
    stop('N.min must be a single value')
  }
  if (!is.null(N.max) && (length(N.max) != 1 || is.na(N.max))) {
    stop('N.max must be a single value')
  }
  check_delta(Delta.A, 'RD', alternative, 'Delta.A', margin, 'RD')
  # The test and the level of the re-estimation are validated as for a rejection region
  ss.Test <- check_paired_args(1, ss.alpha, ss.Test, n.grid, 0, ref.pvalue, alternative,
                               margin)$Test
  # The Berger-Boos parameter matters only for the unconditional test
  if (ss.Test != 'McNemar-uncond') bb.gamma <- 0
  if (restricted && is.null(N)) stop('N must be supplied when restricted is TRUE')
  if (ss.method == 'conditional.power' && is.null(N.max)) {
    stop("ss.method = 'conditional.power' requires N.max")
  }
  lo <- max(n1, if (restricted) N else 0, if (!is.null(N.min)) ceiling(N.min) else 0)
  hi <- if (!is.null(N.max)) floor(N.max) else Inf
  if (hi < lo) {
    stop('N.max is smaller than the smallest final number of pairs the design allows')
  }
  psi.hat <- D1 / n1
  psi <- pmax(psi.hat, abs(Delta.A))
  p10 <- (psi + Delta.A) / 2
  p01 <- (psi - Delta.A) / 2
  n.raw <- rep(NA_real_, length(D1))
  limit <- rep(NA_integer_, length(D1))
  if (ss.method == 'conditional.power') {
    # Probability that a discordant pair favours treatment 1
    pi.dir <- ifelse(psi > 0, p10 / psi, 0.5)
    N.re <- rep(NA_integer_, length(D1))
    todo <- seq_along(D1)
    # Scan upwards from the smallest allowed number of pairs; each region is computed once
    # and serves all numbers of interim discordant pairs not yet settled
    for (n in seq(lo, hi)) {
      rr <- get_paired_rr(n, ss.alpha, ss.Test, alternative, n.grid, bb.gamma, ref.pvalue,
                          margin)
      cp <- paired_cp(rr, n1, as.integer(D1[todo]), psi[todo], pi.dir[todo])
      hit <- cp %>=% tar.power
      N.re[todo[hit]] <- as.integer(n)
      todo <- todo[!hit]
      if (length(todo) == 0) break
    }
    N.re[todo] <- as.integer(hi)
  } else if (ss.method == 'standard') {
    n.raw <- paired_raw_n(p10, p01, ss.alpha, tar.power, alternative, margin)
    N.re <- as.integer(pmax(1, ceil_tol(n.raw)))
  } else {
    key <- paste(p10, p01)
    first <- which(!duplicated(key))
    s <- paired_ss_search(p10[first], p01[first], ss.alpha, tar.power, ss.Test,
                          alternative, n.grid, bb.gamma, ref.pvalue, margin, search,
                          search.limit)
    idx <- match(key, key[first])
    N.re <- s$N[idx]
    limit <- s$limit[idx]
  }
  out <- data.frame(D1 = as.integer(D1), psi.hat = psi.hat, psi.plan = psi,
                    p10.plan = p10, p01.plan = p01, n.raw = n.raw, N.re = N.re,
                    N = as.integer(pmin(hi, pmax(lo, N.re))))
  attr(out, 'N.limit') <- limit
  out
}
