#' Adjusted Significance Level of Blinded Sample Size Re-Estimation for a Stratified Design
#'
#' Lowers the nominal significance level of a two-stratum design with blinded sample size
#' re-estimation until the largest type I error rate over the common response
#' probabilities of the strata does not exceed the target level \code{alpha}, and returns
#' the adjusted level together with that of the fixed design with the planned group sizes.
#'
#' @inheritParams StratTypeIErrorBSSR
#' @param alpha Target level of significance, which the type I error rate must not exceed
#' @param ss.alpha Level of significance used for the re-estimation when
#'   \code{adjust = 'test'}. Default is \code{alpha}. It is ignored when
#'   \code{adjust = 'both'}
#' @param maximize How the largest type I error rate at a level is located, as in
#'   \code{\link{StratTypeIErrorBSSR}}. With \code{'certified'} (default) the adjusted
#'   level controls the type I error rate over the whole rectangle spanned by
#'   \code{theta}. With \code{'grid'} the control is assessed at the points of
#'   \code{theta}
#' @param adjust Which part of the design uses the adjusted level. \code{'test'}
#'   (default) applies it to the final analysis only, so the re-estimation keeps using
#'   \code{ss.alpha}. \code{'both'} applies it to the re-estimation as well
#' @param tol Tolerance of the bisection, relative to \code{alpha}, a single value in
#'   (0, 1). Default is \code{1e-8}
#' @param step Decrement of the level in the search used when \code{adjust = 'both'}.
#'   Default is \code{1e-5}
#'
#' @return An object of class \code{bbssr_alphaadj}, a data frame with one row for the
#'   BSSR design and one for the fixed design containing:
#' \describe{
#'   \item{Design}{\code{'BSSR'} or \code{'TRAD'}}
#'   \item{alpha}{Target level}
#'   \item{max.TIE}{Largest type I error rate at the nominal level \code{alpha}}
#'   \item{alpha.adj}{Adjusted nominal level. It is \code{alpha} when the nominal level
#'     already controls the type I error rate, and 0, with a warning, when no level
#'     examined by the bisection does. The search of \code{adjust = 'both'} for the BSSR
#'     design stops with an error instead when it reaches 0}
#'   \item{max.TIE.adj}{Largest type I error rate at the adjusted level}
#'   \item{theta.S1.adj, theta.S2.adj}{Common response probabilities of the strata at
#'     which \code{max.TIE.adj} occurs, \code{NA} when \code{alpha.adj} is 0}
#'   \item{max.TIE.bound}{Upper bound of the type I error rate at the nominal level with
#'     \code{maximize = 'certified'}, otherwise \code{NA}}
#'   \item{max.TIE.adj.bound}{Upper bound of the type I error rate at the adjusted level
#'     with \code{maximize = 'certified'}, otherwise \code{NA}}
#' }
#'
#' @details
#' The searches are those of \code{\link{BinaryAlphaAdjBSSR}}. With \code{adjust = 'test'}
#' the final group sizes do not depend on the level, the largest type I error rate is a
#' non-decreasing step function of the level, and a bisection locates the level above
#' which it exceeds \code{alpha}, to within \code{tol * alpha}. The reported level is the
#' largest value with at most six significant digits that gives the same critical values.
#' With \code{adjust = 'both'} the level is lowered from \code{alpha} in steps of
#' \code{step} until the type I error rate is controlled. Each level requires the
#' coefficients of the type I error rate of \code{\link{StratTypeIErrorBSSR}}, so the
#' search takes about as long as that many evaluations of the type I error rate.
#'
#' With a random allocation (\code{alloc} other than \code{'fixed'}) the type I error
#' rate of the design with re-estimation at each level is that of
#' \code{\link{StratTypeIErrorBSSR}}, the average over the same \code{nseq} allocation
#' sequences at every level, so each level takes about \code{nseq} times as long as with
#' the fixed allocation. The adjusted level controls this average, whose Monte Carlo error
#' it does not account for. The rate of the fixed design with the planned total and the
#' same allocation is exact. Both levels are reported as the search found them, without
#' the rounding to six significant digits, since the critical values differ between the
#' allocation sequences and between the stratum and group sizes.
#'
#' @examples
#' \donttest{
#' StratAlphaAdjBSSR(Delta.A = 0.25, N1 = c(15, 15), N2 = c(15, 15), omega = 0.5, r = 1,
#'                   alpha = 0.025, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
#'                   theta = seq(0.05, 0.95, by = 0.05))
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratTypeIErrorBSSR}}
#' @export
StratAlphaAdjBSSR <- function(Delta.A, N1, N2, omega = NULL, r, alpha, tar.power, Test,
                              restricted = FALSE,
                              alternative = c('greater', 'less', 'two.sided'),
                              tsmethod = c('minlike', 'central', 'blaker'),
                              effect = c('RD', 'RR', 'OR'),
                              ss.method = c('exact', 'standard'),
                              ss.Test = Test, ss.alpha = alpha,
                              ss.strata = c('stratified', 'pooled'),
                              search = c('crossing', 'smallest', 'stable'),
                              search.limit = c(2, 50), N.min = NULL, N.max = NULL,
                              n.interim = NULL, frac = NULL,
                              theta = seq(0, 1, by = 0.05),
                              maximize = c('certified', 'grid'),
                              adjust = c('test', 'both'), tol = 1e-8, step = 1e-5,
                              alloc = c('fixed', 'stratum', 'group', 'both'), block = 4,
                              ss.frac = c('observed', 'planned'), nseq = 100) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  ss.strata <- match.arg(ss.strata)
  maximize <- match.arg(maximize)
  adjust <- match.arg(adjust)
  alloc <- match.arg(alloc)
  ss.frac <- match.arg(ss.frac)
  if (length(tol) != 1 || !is.numeric(tol) || is.na(tol) || tol <= 0 || tol >= 1) {
    stop('tol must be a single value in (0, 1)')
  }
  if (length(step) != 1 || !is.numeric(step) || is.na(step) || step <= 0 ||
      step >= alpha) {
    stop('step must be a single value in (0, alpha)')
  }
  if (length(N1) != 2 || length(N2) != 2) {
    stop('the adjusted level is available for designs with two strata')
  }
  d <- strat_design(N1, N2, r, frac, omega, n.interim, alpha, Test, alternative, tsmethod,
                    ss.method, 1e5)
  a <- d$a
  pts <- strat_theta(theta, 2)
  box <- apply(pts, 2, range)
  if (maximize == 'certified' && any(box[2, ] - box[1, ] <= 0)) {
    stop("maximize = 'certified' requires theta to span an interval in each stratum")
  }
  re.args <- function(level.ss) {
    list(Delta.A = Delta.A, effect = effect, frac = d$frac, r = d$r,
         ss.method = ss.method, ss.Test = ss.Test, ss.alpha = level.ss,
         tar.power = tar.power, N1.plan = a$N1, N2.plan = a$N2, restricted = restricted,
         N.min = N.min, N.max = N.max, search = search, search.limit = search.limit,
         ss.strata = ss.strata)
  }
  # A design at a level is its type I error rate f at points (x, y), the coefficients of
  # the rate and the final designs. assess() returns the largest rate at the points of
  # theta; with a probe, a point, the rate is first evaluated there, and a value above
  # alpha is returned at once, since the largest rate is at least this value
  make_design <- function(n1, n2, sizes, id, level) {
    designs <- lapply(sizes, function(s) {
      get_strat_crit(s$N1, s$N2, level, a$Test, a$alternative, a$tsmethod)
    })
    coef <- strat_tie_coef(n1, n2, designs, id)
    list(coef = coef, designs = designs,
         f = function(x, y) strat_bern_eval(coef, x, y))
  }
  assess <- function(des, probe = NULL) {
    if (!is.null(probe)) {
      v <- des$f(probe[1], probe[2])
      if (isTRUE(v > alpha)) return(list(x = probe, y = v, bound = NA_real_))
    }
    y <- des$f(pts[, 1], pts[, 2])
    k <- which.max(y)
    list(x = unname(pts[k, ]), y = y[k], bound = NA_real_)
  }
  certify <- if (maximize == 'certified') {
    function(des, m) {
      cm <- strat_tie_certify(des$coef, box)
      if (cm$y > m$y) {
        m$x <- cm$x
        m$y <- cm$y
      }
      m$bound <- max(cm$bound, m$y)
      m
    }
  }
  passes <- function(m) m$y <= alpha && (is.na(m$bound) || m$bound <= alpha)
  # Largest p-value rejected and smallest p-value not rejected at a level, without the
  # values -1 and 2 that stand for none, from which report_level() chooses the level
  brackets <- function(designs, level) {
    lapply(designs, function(s) {
      v <- strat_level_bracket(s$N1, s$N2, a$type, a$alt, a$ts, level)
      v[v >= 0 & v <= 1]
    })
  }
  # With a random allocation the rate of the re-estimation design is the average over the
  # allocation sequences, recomputed at each level, and that of the fixed design the exact
  # sum over the stratum and group sizes; the level of the bisection is reported as it is
  if (alloc != 'fixed') {
    st <- strat_alloc_setup(d, alloc, block, nseq, re.args, ss.frac, restricted, N.min,
                            N.max)
    poly_design <- function(coef) {
      list(coef = coef, f = function(x, y) strat_bern_eval(coef, x, y))
    }
    make_seq <- function(lay, level) {
      poly_design(strat_seq_tie_coef(lay$inputs, a$type, a$alt, a$ts, level, numeric(0),
                                     numeric(0))$coef)
    }
  }
  # Fixed design
  plan <- list(list(N1 = a$N1, N2 = a$N2))
  make_trad <- if (alloc == 'fixed') {
    function(level) make_design(c(0L, 0L), c(0L, 0L), plan, 0L, level)
  } else {
    function(level) {
      poly_design(strat_alloc_tie_coef(st$N.plan, st$al, alloc, level, a$Test,
                                       a$alternative, a$tsmethod)$coef)
    }
  }
  res.trad <- bisect_level(make_trad, assess, certify, alpha, tol)
  # The level 0 of a bisection that found no level is reported as it is
  if (alloc == 'fixed' && res.trad$level > 0 && res.trad$level < alpha) {
    res.trad$level <- report_level(brackets(plan, res.trad$level), res.trad$level)
  }
  # Design with re-estimation
  # Only the final group sizes are needed here, so the critical values of strat_map() are
  # not computed; make_design() computes them at each level
  layout_at <- function(level.ss) {
    sm <- strat_map(d$n1, d$n2, alpha, a$Test, a$alternative, a$tsmethod,
                    re.args(level.ss), crit = FALSE)
    list(sizes = sm$designs, id = sm$id)
  }
  if (adjust == 'test') {
    if (alloc == 'fixed') {
      lay <- layout_at(ss.alpha)
      make_bssr <- function(level) make_design(d$n1, d$n2, lay$sizes, lay$id, level)
    } else {
      lay <- st$layout(ss.alpha)
      make_bssr <- function(level) make_seq(lay, level)
    }
    res.bssr <- bisect_level(make_bssr, assess, certify, alpha, tol)
    if (alloc == 'fixed' && res.bssr$level > 0 && res.bssr$level < alpha) {
      res.bssr$level <- report_level(brackets(lay$sizes, res.bssr$level), res.bssr$level)
    }
  } else {
    des_at <- function(level) {
      if (alloc != 'fixed') return(make_seq(st$layout(level), level))
      lay <- layout_at(level)
      make_design(d$n1, d$n2, lay$sizes, lay$id, level)
    }
    des <- des_at(alpha)
    m0 <- assess(des)
    if (!is.null(certify)) m0 <- certify(des, m0)
    k <- 0L
    level <- alpha
    m <- m0
    while (!passes(m)) {
      # Each level is computed from alpha, so that rounding errors do not accumulate
      k <- k + 1L
      level <- alpha - k * step
      if (level < step / 2) stop('no level above zero controls the type I error rate')
      des <- des_at(level)
      m <- assess(des, m$x)
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
  loc <- function(m) if (length(m$x) == 2) m$x else c(NA_real_, NA_real_)
  out <- data.frame(
    Design = c('BSSR', 'TRAD'),
    alpha = alpha,
    max.TIE = c(res.bssr$m0$y, res.trad$m0$y),
    alpha.adj = c(res.bssr$level, res.trad$level),
    max.TIE.adj = c(res.bssr$m$y, res.trad$m$y),
    theta.S1.adj = c(loc(res.bssr$m)[1], loc(res.trad$m)[1]),
    theta.S2.adj = c(loc(res.bssr$m)[2], loc(res.trad$m)[2]),
    max.TIE.bound = c(res.bssr$m0$bound, res.trad$m0$bound),
    max.TIE.adj.bound = c(res.bssr$m$bound, res.trad$m$bound),
    stringsAsFactors = FALSE
  )
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'adjust') <- adjust
  attr(out, 'maximize') <- maximize
  if (maximize == 'certified') attr(out, 'interval') <- t(box)
  attr(out, 'nuisance') <- 'theta'
  attr(out, 'N1') <- a$N1
  attr(out, 'N2') <- a$N2
  attr(out, 'n1.interim') <- d$n1
  attr(out, 'n2.interim') <- d$n2
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.strata') <- ss.strata
  attr(out, 'search') <- search
  attr(out, 'margin') <- 0
  attr(out, 'margin.scale') <- 'RD'
  attr(out, 'alloc') <- alloc
  if (alloc != 'fixed') {
    attr(out, 'block') <- st$al$block
    attr(out, 'ss.frac') <- ss.frac
    attr(out, 'nseq') <- as.integer(nseq)
  }
  class(out) <- c('bbssr_alphaadj', 'data.frame')
  out
}
