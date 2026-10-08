#' Adjusted Significance Level of a Blinded Sample Size Re-estimation Design
#'
#' Lowers the nominal significance level of a two-arm trial with a binary endpoint and
#' blinded sample size re-estimation (BSSR) until the largest type I error rate over the
#' common response probability, or over the pooled response probability on the boundary
#' of a non-inferiority hypothesis, does not exceed the target level \code{alpha}, and
#' returns the adjusted level together with that of the fixed-sample design.
#'
#' @inheritParams BinaryTypeIErrorBSSR
#' @param alpha Target level of significance, which the type I error rate must not exceed
#' @param ss.alpha Level of significance used for the re-estimation when
#'   \code{adjust = 'test'}. Default is \code{alpha}. It is ignored when
#'   \code{adjust = 'both'}
#' @param adjust Which part of the design uses the adjusted level. \code{'test'}
#'   (default) applies it to the final analysis only, so the re-estimation keeps using
#'   \code{ss.alpha}. \code{'both'} applies it to the re-estimation as well
#' @param tol Tolerance of the bisection, relative to \code{alpha}, a single value in
#'   (0, 1). The bisection is used when \code{adjust = 'test'} and, whatever
#'   \code{adjust}, for the fixed-sample design. Default is \code{1e-8}
#' @param step Decrement of the level in the search used when \code{adjust = 'both'}.
#'   Default is \code{1e-5}
#' @param maximize How the largest type I error rate at a level is located, as in
#'   \code{\link{BinaryTypeIErrorBSSR}}. With \code{'certified'} (default) the adjusted
#'   level controls the type I error rate over the whole interval from the smallest to the
#'   largest value of \code{theta}, restricted with a margin to the values at which both
#'   response probabilities on the null boundary lie in the unit interval, see Details.
#'   With \code{'refined'} and \code{'grid'} the control is assessed on the grid
#'
#' @return An object of class \code{bbssr_alphaadj}, a data frame with one row for the
#'   BSSR design and one for the fixed-sample design containing:
#' \describe{
#'   \item{Design}{\code{'BSSR'} or \code{'TRAD'}}
#'   \item{alpha}{Target level}
#'   \item{max.TIE}{Largest type I error rate at the nominal level \code{alpha}}
#'   \item{alpha.adj}{Adjusted nominal level, see Details. It is \code{alpha} when the
#'     nominal level already controls the type I error rate, and 0, with a warning, when
#'     no level examined does}
#'   \item{max.TIE.adj}{Largest type I error rate at the adjusted level, 0 when
#'     \code{alpha.adj} is 0}
#'   \item{theta.adj}{Common response probability, or pooled response probability on the
#'     null boundary, at which \code{max.TIE.adj} occurs, \code{NA} when
#'     \code{alpha.adj} is 0}
#'   \item{max.TIE.bound}{Upper bound of the type I error rate at the nominal level with
#'     \code{maximize = 'certified'}, otherwise \code{NA}}
#'   \item{max.TIE.adj.bound}{Upper bound of the type I error rate at the adjusted level
#'     with \code{maximize = 'certified'}, otherwise \code{NA}}
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
#' rate is a non-decreasing step function of the level. The bisection locates, to within
#' \code{tol * alpha}, the level above which the largest rate exceeds \code{alpha}. The
#' package rejects a p-value at a level when it is below the level by more than a
#' tolerance of about 1.49e-8, so every level between the largest p-value rejected at the
#' level found plus this tolerance and the smallest p-value not rejected gives the same
#' rejection regions, also when a p-value is rejected if it is below the level or at most
#' equal to it. The reported \code{alpha.adj} is the largest value with at most six
#' significant digits in this range, or with more digits if there is none, so the printed
#' value can be used as it is. The same applies to the fixed-sample design.
#'
#' With \code{adjust = 'both'} the re-estimated sample sizes change with the level as
#' well, so the type I error rate need not be monotone in the level. The level is then
#' lowered from \code{alpha} in steps of \code{step} until the type I error rate is
#' controlled. The result is the first level among \code{alpha - step},
#' \code{alpha - 2 * step}, and so on, at which the type I error rate is controlled, and
#' it can change with \code{step}. The levels that control the type I error rate need not
#' form an interval, so a level below the result can fail, and a level above it that is
#' not on this grid can pass. This search re-estimates the sample size at every step and
#' can take much longer than the bisection. The fixed-sample design always uses the
#' bisection, since its sample size does not depend on the level. The script
#' \code{reproduce-published.R} in the folder
#' \code{system.file('reproduce', package = 'bbssr')} reproduces the published values of
#' Section 5 of Friede and Kieser (2004) and of Example 21.1 of Kieser (2020) only when the
#' re-estimation also uses the adjusted level, as with \code{adjust = 'both'}.
#'
#' The type I error rate can also be controlled within a confidence interval for the
#' common response probability, or for the pooled response probability on the null
#' boundary, instead of over its whole range (Kieser, 2020, Section 22.1.3). A value
#' \code{gamma} between 0 and \code{alpha} is fixed in advance, and the final analysis
#' uses a level at which the largest type I error rate over the \code{1 - gamma}
#' confidence interval computed from the data of the completed trial does not exceed
#' \code{alpha - gamma}. The type I error rate is then controlled at \code{alpha}. This
#' level is found when \code{theta} spans the interval and \code{alpha} is set to
#' \code{alpha - gamma}. The interval is known only at the end of the trial, so the
#' re-estimation keeps the nominal level, with \code{adjust = 'test'} and \code{ss.alpha}
#' set to the nominal level. Neither search examines levels above \code{alpha - gamma}.
#' The script \code{reproduce-published.R} reproduces the adjusted level of Example 23.1
#' of Kieser (2020), which uses the Clopper-Pearson interval with \code{gamma = 1e-4}.
#'
#' With \code{maximize = 'certified'} the levels are assessed on the grid with the
#' refinement of \code{maximize = 'refined'}, and the level found is then certified as in
#' \code{\link{BinaryTypeIErrorBSSR}}: it is accepted only if the upper bound of the type I
#' error rate over the interval from the smallest to the largest value of \code{theta}
#' does not exceed \code{alpha}. If the bound exceeds \code{alpha}, the grid has missed a
#' value of \code{theta} at which the level fails, or the largest rate lies within 1e-12
#' below \code{alpha}, and the search continues below that level. The bisection then
#' certifies every level that passes the assessment on the grid, and the search of
#' \code{adjust = 'both'} certifies every such level from the start. The adjusted level
#' therefore controls the type I error rate over the whole interval. If a certification
#' stops at its limit of 1e5 subdivisions, with a warning, the bound can exceed the
#' largest rate by more than 1e-12, and the level found can be lower than necessary.
#'
#' In both searches the type I error rate at a new level is first evaluated at the single
#' value of \code{theta} where the largest rate of the last level that failed was found.
#' If it exceeds \code{alpha} there, the level fails without the evaluation over the whole
#' grid. The largest rate over \code{theta} is never below the rate at any one value, so
#' such a level also fails on the largest rate. With \code{maximize = 'grid'} this value
#' is a grid point, and every decision is that of the evaluation on the grid. With
#' \code{'refined'} and \code{'certified'} it can lie between the grid points, where the
#' grid can miss a rate above \code{alpha}, so a level can fail that the evaluation on the
#' grid would accept. This can only lower the level found, and with \code{'certified'}
#' the level found is certified in either case.
#'
#' @references
#' Kieser M, Friede T (2000). Re-calculating the sample size in internal pilot study
#' designs with control of the type I error rate. \emph{Statistics in Medicine}, 19(7),
#' 901-911.
#'
#' Friede T, Kieser M (2004). Sample size recalculation for binary data in internal pilot
#' study designs. \emph{Pharmaceutical Statistics}, 3(4), 269-279.
#'
#' Kieser M (2020). \emph{Methods and Applications of Sample Size Calculation and
#' Recalculation in Clinical Trials}. Springer, Cham.
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
                               tsmethod = c('minlike', 'central', 'blaker'),
                               n.grid = 100, bb.gamma = 0,
                               effect = c('RD', 'RR', 'OR'),
                               ss.method = c('exact', 'standard', 'null.variance',
                                             'alternative.variance'),
                               ss.Test = Test, ss.alpha = alpha,
                               rounding = c('group', 'friede-kieser', 'total', 'nearest'),
                               search = c('crossing', 'smallest', 'stable'),
                               search.limit = c(2, 50),
                               N.min = NULL, N.max = NULL, n.interim = NULL,
                               theta = seq(0, 1, by = 0.005),
                               maximize = c('certified', 'refined', 'grid'),
                               adjust = c('test', 'both'), tol = 1e-8, step = 1e-5,
                               margin = 0, ref.pvalue = FALSE,
                               margin.scale = c('RD', 'RR')) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  search <- match.arg(search)
  adjust <- match.arg(adjust)
  maximize <- match.arg(maximize)
  margin.scale <- match.arg(margin.scale)
  if (length(theta) < 1 || anyNA(theta) || any(theta < 0 | theta > 1)) {
    stop('theta must be a vector of values in [0, 1]')
  }
  if (length(tol) != 1 || is.na(tol) || tol <= 0 || tol >= 1) {
    stop('tol must be a single value in (0, 1)')
  }
  if (length(step) != 1 || is.na(step) || step <= 0 || step >= alpha) {
    stop('step must be a single value in (0, alpha)')
  }
  theta <- sort(unique(theta))
  Test <- check_rr_args(N1, N2, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                        margin, margin.scale)$Test
  # Response probabilities on the boundary of the null hypothesis, which are both theta
  # when the margin is 0
  ok <- null_boundary(theta, r, alternative, margin, margin.scale)$ok
  if (!any(ok)) {
    stop('no value of theta gives response probabilities in [0, 1] on the null boundary')
  }
  # Interval over which the largest rate is certified
  interval <- null_range(theta, r, alternative, margin, margin.scale)
  theta <- theta[ok]
  # p-values of a rejection region, from which the region at any level follows
  pvalues <- function(n1, n2) {
    a <- check_rr_args(n1, n2, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                       margin, margin.scale)
    get_pvalue(a$N1, a$N2, a$Test, alternative, tsmethod, a$n.grid, bb.gamma,
               a$ref.pvalue, a$margin, a$margin.scale)
  }
  # A design at a given level is a list with its type I error rate f as a function of
  # theta, its interim setup and its rejection regions. max_tie returns the largest rate
  # over the grid, refined between the grid points unless maximize = 'grid'. With a probe,
  # a value of theta, the rate is first evaluated there, and a value above alpha is
  # returned at once: the largest rate is at least this value, so the level fails either
  # way
  max_tie <- function(des, probe = NULL) {
    if (!is.null(probe)) {
      v <- des$f(probe)
      if (isTRUE(v > alpha)) return(list(x = probe, y = v, bound = NA_real_))
    }
    y <- des$f(theta)
    m <- if (maximize == 'grid') {
      list(x = theta[which.max(y)], y = max(y))
    } else {
      refine_max(des$f, theta, y)
    }
    m$bound <- NA_real_
    m
  }
  # Certified maximum over the interval, combined with the result m of max_tie
  certify <- if (maximize == 'certified') {
    function(des, m) {
      cm <- tie_certify(tie_weights(des$setup, des$rr.list), interval, r, alternative,
                        margin, margin.scale = margin.scale)
      if (cm$y > m$y) {
        m$x <- cm$x
        m$y <- cm$y
      }
      m$bound <- max(cm$bound, m$y)
      m
    }
  }
  passes <- function(m) m$y <= alpha && (is.na(m$bound) || m$bound <= alpha)
  # Fixed-sample design, a design without an interim stage
  pv.fixed <- pvalues(N1, N2)
  setup.fixed <- fixed_setup(N1, N2)
  make_trad <- function(level) {
    rr <- pv.fixed %<<% level
    f <- function(t) {
      b <- null_boundary(t, r, alternative, margin, margin.scale)
      vapply(seq_along(t), function(i) {
        power_from_rr(rr, dbinom(0:N1, N1, b$p1[i]), dbinom(0:N2, N2, b$p2[i]))
      }, numeric(1))
    }
    list(f = f, setup = setup.fixed, rr.list = list(rr))
  }
  res.trad <- bisect_level(make_trad, max_tie, certify, alpha, tol)
  if (res.trad$level < alpha) res.trad$level <- report_level(list(pv.fixed), res.trad$level)
  # BSSR design
  design_at <- function(level.ss) {
    map <- bssr_map(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                    restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                    ss.method, ss.Test, level.ss, rounding, N.min, N.max, ref.pvalue,
                    margin, search = search, search.limit = search.limit,
                    margin.scale = margin.scale)
    setup <- bssr_setup(map)
    pv.list <- lapply(seq_along(setup$N1), function(k) pvalues(setup$N1[k], setup$N2[k]))
    list(setup = setup, pv.list = pv.list)
  }
  make_bssr <- function(d, level) {
    rr.list <- lapply(d$pv.list, function(m) m %<<% level)
    f <- function(t) {
      b <- null_boundary(t, r, alternative, margin, margin.scale)
      bssr_reject(d$setup, rr.list, b$p1, b$p2)
    }
    list(f = f, setup = d$setup, rr.list = rr.list)
  }
  if (adjust == 'test') {
    d <- design_at(ss.alpha)
    res.bssr <- bisect_level(function(level) make_bssr(d, level), max_tie, certify, alpha,
                             tol)
    if (res.bssr$level < alpha) res.bssr$level <- report_level(d$pv.list, res.bssr$level)
  } else {
    # The nominal level is certified, and so is every lower level that passes the
    # assessment on the grid
    des <- make_bssr(design_at(alpha), alpha)
    m0 <- max_tie(des)
    if (!is.null(certify)) m0 <- certify(des, m0)
    k <- 0L
    level <- alpha
    m <- m0
    while (!passes(m)) {
      # Each level is computed from alpha, so that rounding errors do not accumulate
      k <- k + 1L
      level <- alpha - k * step
      if (level < step / 2) stop('no level above zero controls the type I error rate')
      des <- make_bssr(design_at(level), level)
      m <- max_tie(des, m$x)
      if (!is.null(certify) && m$y <= alpha) m <- certify(des, m)
    }
    res.bssr <- list(level = level, m0 = m0, m = m)
  }
  zero <- c(BSSR = res.bssr$level, TRAD = res.trad$level) == 0
  if (any(zero)) {
    warning('no level examined controls the type I error rate of the ',
            paste(names(zero)[zero], collapse = ' and '),
            ' design, and the adjusted level is reported as 0')
  }
  out <- data.frame(
    Design = c('BSSR', 'TRAD'),
    alpha = alpha,
    max.TIE = c(res.bssr$m0$y, res.trad$m0$y),
    alpha.adj = c(res.bssr$level, res.trad$level),
    max.TIE.adj = c(res.bssr$m$y, res.trad$m$y),
    theta.adj = c(res.bssr$m$x, res.trad$m$x),
    max.TIE.bound = c(res.bssr$m0$bound, res.trad$m0$bound),
    max.TIE.adj.bound = c(res.bssr$m$bound, res.trad$m$bound),
    stringsAsFactors = FALSE
  )
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'adjust') <- adjust
  attr(out, 'maximize') <- maximize
  if (maximize == 'certified') attr(out, 'interval') <- interval
  attr(out, 'N1') <- N1
  attr(out, 'N2') <- N2
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'search') <- search
  attr(out, 'margin') <- margin
  attr(out, 'margin.scale') <- margin.scale
  attr(out, 'ref.pvalue') <- ref.pvalue
  class(out) <- c('bbssr_alphaadj', 'data.frame')
  out
}
