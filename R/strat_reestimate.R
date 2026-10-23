#' Re-Estimated Group Sizes of a Stratified Design from Blinded Interim Data
#'
#' Internal helper returning the final group sizes of the strata of a stratified design
#' with blinded sample size re-estimation, for each row of the interim pooled numbers of
#' responders \code{s}. In stratum \code{j} the blinded estimate of the pooled response
#' probability is \code{s[, j] / (n1[j] + n2[j])}, which \code{split_pooled} divides
#' between the groups with the assumed effect \code{Delta.A}, common to the strata on the
#' scale \code{effect}, and the allocation ratio of the stratum.
#'
#' With \code{ss.method = 'standard'} the nominal total sample size is that of
#' \code{strat_raw_n}, rounded up, with the probabilities before truncation to the unit
#' interval, which keep an assumed risk difference. With \code{ss.method = 'exact'} it is
#' the result of \code{strat_ss_search} for the truncated probabilities. The final group
#' sizes are those of \code{strat_sizes} for the nominal size, kept within the range of
#' \code{strat_bounds}: at least the interim sizes, at least the planned sizes when
#' \code{restricted = TRUE}, a total of at least \code{N.min} and at most \code{N.max}.
#' A row whose probabilities show no effect in the direction of the alternative admits no
#' sample size, and keeps the planned group sizes, raised to the interim sizes where
#' these are larger.
#'
#' @param s Interim pooled numbers of responders, a matrix with one row per interim
#'   outcome and one column per stratum
#' @param n1 Interim sizes of group 1 in the strata
#' @param n2 Interim sizes of group 2 in the strata
#' @param Delta.A Assumed effect, common to the strata
#' @param effect \code{'RD'}, \code{'RR'} or \code{'OR'}
#' @param frac Fractions of the strata, summing to 1
#' @param r Allocation ratios of the strata
#' @param ss.method \code{'exact'} or \code{'standard'}
#' @param ss.Test Test of the re-estimation
#' @param ss.alpha Level of the re-estimation
#' @param tar.power Target power
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod Two-sided method of the stratified exact test
#' @param N1.plan Planned sizes of group 1 in the strata, or \code{NULL}
#' @param N2.plan Planned sizes of group 2 in the strata, or \code{NULL}
#' @param restricted Logical. Whether the final sizes are at least the planned sizes
#' @param N.min Smallest final total sample size, or \code{NULL}
#' @param N.max Largest final total sample size, or \code{NULL}
#' @param search Search of \code{strat_ss_search}
#' @param search.limit Limit of the search \code{'stable'}
#'
#' @return A list with the matrices \code{hat.p}, \code{hat.p1} and \code{hat.p2} of the
#'   blinded estimates and of the truncated planning probabilities, the vectors
#'   \code{n.raw} (normal approximation, \code{NA} under \code{ss.method = 'exact'}),
#'   \code{n.re} (nominal size before the bounds, \code{NA} for a row without a sample
#'   size), \code{limit} (limit of the search \code{'stable'}, otherwise \code{NA}), and
#'   the integer matrices \code{N1} and \code{N2} of the final group sizes
#'
#' @keywords internal
#' @noRd
strat_reestimate <- function(s, n1, n2, Delta.A, effect, frac, r, ss.method, ss.Test,
                             ss.alpha, tar.power, alternative, tsmethod, N1.plan, N2.plan,
                             restricted, N.min, N.max, search, search.limit) {
  J <- length(n1)
  s <- matrix(as.double(s), ncol = J)
  check_delta(Delta.A, effect, alternative, 'Delta.A', 0, 'RD')
  if (length(tar.power) != 1 || is.na(tar.power) || tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  if (restricted && (is.null(N1.plan) || is.null(N2.plan))) {
    stop('the planned sizes N1 and N2 must be supplied when restricted is TRUE')
  }
  a <- check_strat_args(rep(1, J), rep(1, J), ss.alpha, ss.Test, alternative, tsmethod)
  hat.p <- sweep(s, 2, n1 + n2, '/')
  sp1 <- sp2 <- matrix(NA_real_, nrow(s), J)
  for (j in seq_len(J)) {
    sp <- split_pooled(hat.p[, j], Delta.A, r[j], effect)
    sp1[, j] <- sp$p1
    sp2[, j] <- sp$p2
  }
  # The matrix comes first in pmin() and pmax(), which keep its dimensions
  hat.p1 <- pmin(pmax(sp1, 0), 1)
  hat.p2 <- pmin(pmax(sp2, 0), 1)
  n.raw <- rep(NA_real_, nrow(s))
  limit <- rep(NA_integer_, nrow(s))
  n.re <- rep(NA_real_, nrow(s))
  if (ss.method == 'standard') {
    n.raw <- strat_raw_n(sp1, sp2, frac, r, ss.alpha, tar.power, a$alternative, a$cc)
    ok <- is.finite(n.raw)
    n.re[ok] <- pmax(1, ceil_tol(n.raw[ok]))
  } else {
    ok <- is.finite(strat_raw_n(hat.p1, hat.p2, frac, r, ss.alpha, tar.power,
                                a$alternative, a$cc))
    key <- apply(cbind(hat.p1, hat.p2), 1, paste, collapse = '_')
    first <- which(!duplicated(key) & ok)
    if (length(first) > 0) {
      res <- strat_ss_search(hat.p1[first, , drop = FALSE], hat.p2[first, , drop = FALSE],
                             frac, r, ss.alpha, tar.power, a$Test, a$alternative,
                             a$tsmethod, search, search.limit)
      idx <- match(key, key[first])
      n.re[ok] <- res$n[idx[ok]]
      limit[ok] <- res$limit[idx[ok]]
    }
  }
  if (any(!ok) && (is.null(N1.plan) || is.null(N2.plan))) {
    stop('the recovered response probabilities show no effect in the direction of the ',
         'alternative, so no sample size can be re-estimated; supply the planned sizes ',
         'N1 and N2')
  }
  low1 <- if (restricted) pmax(n1, N1.plan) else n1
  low2 <- if (restricted) pmax(n2, N2.plan) else n2
  b <- strat_bounds(frac, r, low1, low2, N.min, N.max)
  n.fin <- pmin(b[2], pmax(b[1], n.re))
  N1 <- matrix(NA_integer_, nrow(s), J)
  N2 <- matrix(NA_integer_, nrow(s), J)
  if (any(ok)) {
    sz <- strat_sizes(n.fin[ok], frac, r)
    N1[ok, ] <- sz$N1
    N2[ok, ] <- sz$N2
  }
  if (any(!ok)) {
    N1[!ok, ] <- rep(as.integer(pmax(N1.plan, n1)), each = sum(!ok))
    N2[!ok, ] <- rep(as.integer(pmax(N2.plan, n2)), each = sum(!ok))
  }
  list(hat.p = hat.p, hat.p1 = hat.p1, hat.p2 = hat.p2, n.raw = n.raw,
       n.re = n.re, limit = limit, N1 = N1, N2 = N2)
}
