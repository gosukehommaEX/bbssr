#' Adjusted Significance Level of a Blinded Sample Size Re-estimation Design
#'
#' Finds the largest nominal significance level at which a two-arm trial with a binary
#' endpoint and blinded sample size re-estimation (BSSR) keeps its type I error rate at or
#' below the target level \code{alpha} for every common response probability, together
#' with the corresponding level of the fixed-sample design.
#'
#' @inheritParams BinaryTypeIErrorBSSR
#' @param alpha Target level of significance, which the type I error rate must not exceed
#' @param ss.alpha Level of significance used for the re-estimation when
#'   \code{adjust = 'test'}. Default is \code{alpha}. It is ignored when
#'   \code{adjust = 'both'}
#' @param adjust Which part of the design uses the adjusted level. \code{'test'}
#'   (default) applies it to the final analysis only, so the re-estimation keeps using
#'   \code{ss.alpha}. \code{'both'} applies it to the re-estimation as well
#' @param tol Tolerance of the bisection used when \code{adjust = 'test'}, relative to
#'   \code{alpha}. Default is \code{1e-8}
#' @param step Decrement of the level in the search used when \code{adjust = 'both'}.
#'   Default is \code{1e-5}
#'
#' @return An object of class \code{bbssr_alphaadj}, a data frame with one row for the
#'   BSSR design and one for the fixed-sample design containing:
#' \describe{
#'   \item{Design}{\code{'BSSR'} or \code{'TRAD'}}
#'   \item{alpha}{Target level}
#'   \item{max.TIE}{Largest type I error rate at the nominal level \code{alpha}}
#'   \item{alpha.adj}{Adjusted nominal level}
#'   \item{max.TIE.adj}{Largest type I error rate at the adjusted level}
#'   \item{theta.adj}{Common response probability at which \code{max.TIE.adj} occurs}
#' }
#'
#' @details
#' Following Kieser and Friede (2000) and Friede and Kieser (2004), the nominal level of
#' the final analysis is lowered until the largest type I error rate over the common
#' response probability, evaluated as in \code{\link{BinaryTypeIErrorBSSR}}, does not
#' exceed \code{alpha}.
#'
#' With \code{adjust = 'test'} the final sample sizes do not depend on the nominal level,
#' and the rejection regions shrink as the level decreases, so the largest type I error
#' rate is a non-decreasing step function of the level. The adjusted level is then found
#' by bisection to within \code{tol * alpha}, and the reported value is the largest level
#' examined that controls the type I error rate.
#'
#' With \code{adjust = 'both'} the re-estimated sample sizes change with the level as
#' well, so the type I error rate need not be monotone in the level. The level is then
#' lowered from \code{alpha} in steps of \code{step} until the type I error rate is
#' controlled. This search re-estimates the sample
#' size at every step and can take much longer than the bisection. The fixed-sample design
#' always uses the bisection, since its sample size does not depend on the level.
#'
#' @references
#' Kieser M, Friede T (2000). Re-calculating the sample size in internal pilot study
#' designs with control of the type I error rate. \emph{Statistics in Medicine}, 19(7),
#' 901-911.
#'
#' Friede T, Kieser M (2004). Sample size recalculation for binary data in internal pilot
#' study designs. \emph{Pharmaceutical Statistics}, 3(4), 269-279.
#'
#' @examples
#' \donttest{
#' BinaryAlphaAdjBSSR(
#'   Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
#'   theta = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryTypeIErrorBSSR}}
#' @export
#' @import fpCompare
#' @importFrom stats dbinom
BinaryAlphaAdjBSSR <- function(Delta.A, N1, N2, omega = NULL, r, alpha, tar.power, Test,
                               restricted = FALSE,
                               alternative = c('greater', 'less', 'two.sided'),
                               tsmethod = c('minlike', 'central'),
                               n.grid = 100, bb.gamma = 0,
                               effect = c('RD', 'RR', 'OR'),
                               ss.method = c('exact', 'standard', 'null.variance'),
                               ss.Test = Test, ss.alpha = alpha,
                               rounding = c('group', 'friede-kieser', 'total'),
                               N.min = NULL, N.max = NULL, n.interim = NULL,
                               theta = seq(0.005, 0.995, by = 0.005),
                               adjust = c('test', 'both'), tol = 1e-8, step = 1e-5,
                               ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  adjust <- match.arg(adjust)
  if (length(theta) < 1 || anyNA(theta) || any(theta < 0 | theta > 1)) {
    stop('theta must be a vector of values in [0, 1]')
  }
  if (length(tol) != 1 || is.na(tol) || tol <= 0) stop('tol must be a single positive value')
  if (length(step) != 1 || is.na(step) || step <= 0 || step >= alpha) {
    stop('step must be a single value in (0, alpha)')
  }
  theta <- sort(unique(theta))
  Test <- check_rr_args(N1, N2, alpha, Test, n.grid, bb.gamma, ref.pvalue)$Test
  # p-values of a rejection region, from which the region at any level follows
  pvalues <- function(n1, n2) {
    a <- check_rr_args(n1, n2, alpha, Test, n.grid, bb.gamma, ref.pvalue)
    get_pvalue(a$N1, a$N2, a$Test, alternative, tsmethod, a$n.grid, bb.gamma,
               a$ref.pvalue)
  }
  # Largest type I error rate over theta of a design whose rejection regions are given
  max_tie <- function(f) refine_max(f, theta, f(theta))
  # Bisection over the level for designs whose sample sizes do not depend on the level
  bisect <- function(make_f) {
    m <- max_tie(make_f(alpha))
    m0 <- m
    if (m$y <= alpha) return(list(level = alpha, m0 = m0, m = m))
    lo <- 0
    hi <- alpha
    m.lo <- list(x = NA_real_, y = 0)
    while (hi - lo > tol * alpha) {
      mid <- (lo + hi) / 2
      m <- max_tie(make_f(mid))
      if (m$y <= alpha) {
        lo <- mid
        m.lo <- m
      } else {
        hi <- mid
      }
    }
    list(level = lo, m0 = m0, m = m.lo)
  }
  # Fixed-sample design
  pv.fixed <- pvalues(N1, N2)
  make_trad <- function(level) {
    rr <- pv.fixed %<<% level
    function(t) {
      vapply(t, function(u) power_from_rr(rr, dbinom(0:N1, N1, u), dbinom(0:N2, N2, u)),
             numeric(1))
    }
  }
  res.trad <- bisect(make_trad)
  # BSSR design
  design_at <- function(level.ss) {
    map <- bssr_map(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                    restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                    ss.method, ss.Test, level.ss, rounding, N.min, N.max, ref.pvalue)
    setup <- bssr_setup(map)
    pv.list <- lapply(seq_along(setup$N1), function(k) pvalues(setup$N1[k], setup$N2[k]))
    list(setup = setup, pv.list = pv.list)
  }
  make_bssr <- function(d, level) {
    rr.list <- lapply(d$pv.list, function(m) m %<<% level)
    function(t) bssr_reject(d$setup, rr.list, t, t)
  }
  if (adjust == 'test') {
    d <- design_at(ss.alpha)
    res.bssr <- bisect(function(level) make_bssr(d, level))
  } else {
    level <- alpha
    m0 <- max_tie(make_bssr(design_at(alpha), alpha))
    m <- m0
    while (m$y > alpha) {
      level <- level - step
      if (level <= 0) stop('no level above zero controls the type I error rate')
      m <- max_tie(make_bssr(design_at(level), level))
    }
    res.bssr <- list(level = level, m0 = m0, m = m)
  }
  out <- data.frame(
    Design = c('BSSR', 'TRAD'),
    alpha = alpha,
    max.TIE = c(res.bssr$m0$y, res.trad$m0$y),
    alpha.adj = c(res.bssr$level, res.trad$level),
    max.TIE.adj = c(res.bssr$m$y, res.trad$m$y),
    theta.adj = c(res.bssr$m$x, res.trad$m$x),
    stringsAsFactors = FALSE
  )
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'adjust') <- adjust
  attr(out, 'N1') <- N1
  attr(out, 'N2') <- N2
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ref.pvalue') <- ref.pvalue
  class(out) <- c('bbssr_alphaadj', 'data.frame')
  out
}
