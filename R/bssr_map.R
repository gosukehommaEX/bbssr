#' Final Sample Size for Every Interim Outcome of a Re-Estimation Design
#'
#' Internal helper shared by the functions that evaluate a design with blinded sample
#' size re-estimation. The re-estimated sample size depends on the interim data only
#' through the pooled number of responders \code{s}, so the design is described by one row
#' per value of \code{s}.
#'
#' @inheritParams BinaryPowerBSSR
#'
#' @return A data frame with columns \code{s}, \code{hat.p}, \code{hat.p1},
#'   \code{hat.p2}, \code{n.raw}, \code{N1.re}, \code{N2.re}, \code{N1}, \code{N2} and
#'   \code{N2.limit}, where \code{N1} and \code{N2} are the final sample sizes and
#'   \code{N2.limit} is the limit of \code{search = 'stable'} (otherwise \code{NA}), and
#'   with the interim sample sizes as the attributes \code{n11} and \code{n12}
#'
#' @keywords internal
#' @noRd
bssr_map <- function(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                     restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                     ss.method, ss.Test, ss.alpha, rounding, N.min, N.max, ref.pvalue,
                     margin, search = 'crossing', search.limit = c(2, 50), margin.scale) {
  if (!is.null(n.interim)) {
    if (!is.null(omega)) stop('supply either omega or n.interim, not both')
    if (length(n.interim) != 2 || anyNA(n.interim) || any(n.interim < 1) ||
        any(n.interim != round(n.interim))) {
      stop('n.interim must be two positive integers, the interim sizes of group 1 and 2')
    }
    N11 <- as.integer(n.interim[1])
    N12 <- as.integer(n.interim[2])
  } else {
    if (is.null(omega)) stop('supply either omega or n.interim')
    if (length(omega) != 1 || is.na(omega) || omega <= 0 || omega > 1) {
      stop('omega must be a single value in (0, 1]')
    }
    # The experimental group is obtained from the control group by a single rounding, so
    # the interim allocation is as close to r to 1 as integers allow
    N12 <- as.integer(ceiling(omega * N2))
    N11 <- as.integer(ceiling(r * N12))
  }
  if (length(r) != 1 || is.na(r) || r <= 0) stop('r must be a single positive value')
  if (length(tar.power) != 1 || is.na(tar.power) || tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  if (ss.method == 'exact' && rounding != 'group') {
    stop("rounding must be 'group' when ss.method is 'exact'")
  }
  if (!is.null(N.min) && (length(N.min) != 1 || is.na(N.min))) {
    stop('N.min must be a single value')
  }
  if (!is.null(N.max) && (length(N.max) != 1 || is.na(N.max))) {
    stop('N.max must be a single value')
  }
  check_delta(Delta.A, effect, alternative, 'Delta.A', margin, margin.scale)
  n1 <- N11 + N12
  s <- 0:n1
  hat.p <- s / n1
  sp <- split_pooled(hat.p, Delta.A, r, effect)
  hat.p1 <- pmin(1, pmax(0, sp$p1))
  hat.p2 <- pmin(1, pmax(0, sp$p2))
  # The exact search needs probabilities inside the unit interval. The normal
  # approximation uses the untruncated probabilities, which keep the assumed effect, as in
  # formula (2) of Friede and Kieser (2004)
  ss.p <- if (ss.method == 'exact') list(p1 = hat.p1, p2 = hat.p2) else sp
  re <- reestimate(ss.p$p1, ss.p$p2, r, ss.alpha, tar.power, ss.Test, alternative,
                   tsmethod, n.grid, bb.gamma, ss.method, rounding, N1, N2, ref.pvalue,
                   margin, search, search.limit, margin.scale)
  fin <- final_sizes(re, N11, N12, r, rounding, restricted, N1, N2, N.min, N.max)
  out <- data.frame(s = s, hat.p = hat.p, hat.p1 = hat.p1, hat.p2 = hat.p2,
                    n.raw = re$n.raw, N1.re = re$N1.re, N2.re = re$N2.re,
                    N1 = fin$N1, N2 = fin$N2, N2.limit = re$N2.limit)
  attr(out, 'n11') <- N11
  attr(out, 'n12') <- N12
  out
}
