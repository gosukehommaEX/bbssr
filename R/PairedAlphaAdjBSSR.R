#' Adjusted Significance Level of Blinded Sample Size Re-Estimation for Paired Binary Data
#'
#' Lowers the nominal significance level of a matched-pairs design with blinded sample
#' size re-estimation until the largest type I error rate over the probability of a
#' discordant pair on the boundary of the null hypothesis does not exceed the target level
#' \code{alpha}, and returns the adjusted level together with that of the fixed design.
#'
#' @inheritParams PairedTypeIErrorBSSR
#' @param alpha Target level of significance, which the type I error rate must not exceed
#' @param ss.alpha Level of significance used for the re-estimation when
#'   \code{adjust = 'test'}. Default is \code{alpha}. It is ignored when
#'   \code{adjust = 'both'}
#' @param adjust Which part of the design uses the adjusted level. \code{'test'}
#'   (default) applies it to the final analysis only, so the re-estimation keeps using
#'   \code{ss.alpha}. \code{'both'} applies it to the re-estimation as well
#' @param tol Tolerance of the bisection, relative to \code{alpha}, a single value in
#'   (0, 1). Default is \code{1e-8}
#' @param step Decrement of the level in the search used when \code{adjust = 'both'}.
#'   Default is \code{1e-5}
#' @param maximize How the largest type I error rate at a level is located, as in
#'   \code{\link{PairedTypeIErrorBSSR}}. With \code{'certified'} (default) the adjusted
#'   level controls the type I error rate over the whole interval spanned by the values of
#'   \code{psi} on the boundary. With \code{'refined'} and \code{'grid'} the control is
#'   assessed on the grid
#'
#' @return An object of class \code{bbssr_alphaadj}, a data frame with one row for the
#'   BSSR design and one for the fixed design containing:
#' \describe{
#'   \item{Design}{\code{'BSSR'} or \code{'TRAD'}}
#'   \item{alpha}{Target level}
#'   \item{max.TIE}{Largest type I error rate at the nominal level \code{alpha}}
#'   \item{alpha.adj}{Adjusted nominal level. It is \code{alpha} when the nominal level
#'     already controls the type I error rate, and 0, with a warning, when no level
#'     examined does}
#'   \item{max.TIE.adj}{Largest type I error rate at the adjusted level}
#'   \item{psi.adj}{Probability of a discordant pair at which \code{max.TIE.adj} occurs,
#'     \code{NA} when \code{alpha.adj} is 0}
#'   \item{max.TIE.bound}{Upper bound of the type I error rate at the nominal level with
#'     \code{maximize = 'certified'}, otherwise \code{NA}}
#'   \item{max.TIE.adj.bound}{Upper bound of the type I error rate at the adjusted level
#'     with \code{maximize = 'certified'}, otherwise \code{NA}}
#' }
#' For \code{alternative = 'equivalence'} the column \code{boundary.adj} gives the
#' boundary, \code{'lower'} or \code{'upper'}, on which \code{max.TIE.adj} occurs.
#'
#' @details
#' The searches are those of \code{\link{BinaryAlphaAdjBSSR}}. With
#' \code{adjust = 'test'} the final numbers of pairs do not depend on the nominal level of
#' the final analysis, the largest type I error rate is a non-decreasing function of the
#' level, and the level is found by bisection to within \code{tol * alpha}; the reported
#' level is the largest value with at most six significant digits that gives the same
#' rejection regions. With \code{adjust = 'both'} the level is lowered from \code{alpha}
#' in steps of \code{step} until the type I error rate is controlled. The fixed design
#' always uses the bisection. With \code{maximize = 'certified'} a level is accepted only
#' when the upper bound of the type I error rate over the interval, computed as in
#' \code{\link{PairedTypeIErrorBSSR}}, does not exceed \code{alpha}.
#'
#' The exact conditional McNemar test keeps its level whatever blinded rule is used
#' (Schepers et al., 2026), so it needs no adjustment.
#'
#' @references
#' Kieser M, Friede T (2000). Re-calculating the sample size in internal pilot study
#' designs with control of the type I error rate. \emph{Statistics in Medicine}, 19(7),
#' 901-911.
#'
#' Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026). Blinded sample size
#' review for McNemar's test based on primary and surrogate endpoints. arXiv:2608.17784.
#'
#' @examples
#' \donttest{
#' PairedAlphaAdjBSSR(Delta.A = 0.2, N = 40, n.interim = 20, alpha = 0.05,
#'                    tar.power = 0.8, Test = 'McNemar', alternative = 'two.sided',
#'                    ss.method = 'standard', psi = seq(0.05, 0.95, by = 0.05))
#' }
#'
#' @export
PairedAlphaAdjBSSR <- function(Delta.A, N, omega = NULL, alpha, tar.power, Test,
                               restricted = FALSE,
                               alternative = c('greater', 'less', 'two.sided',
                                               'equivalence'),
                               n.grid = 100, bb.gamma = 0,
                               ss.method = c('exact', 'standard', 'conditional.power'),
                               ss.Test = Test, ss.alpha = alpha,
                               search = c('crossing', 'smallest', 'stable'),
                               search.limit = c(2, 50), N.min = NULL, N.max = NULL,
                               n.interim = NULL, psi = seq(0, 1, by = 0.01),
                               maximize = c('certified', 'refined', 'grid'),
                               adjust = c('test', 'both'), tol = 1e-8, step = 1e-5,
                               margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  maximize <- match.arg(maximize)
  adjust <- match.arg(adjust)
  margin <- norm_margin(alternative, margin, 'RD')
  if (length(psi) < 1 || !is.numeric(psi) || anyNA(psi) || any(psi < 0 | psi > 1)) {
    stop('psi must be a vector of values in [0, 1]')
  }
  if (length(tol) != 1 || is.na(tol) || tol <= 0 || tol >= 1) {
    stop('tol must be a single value in (0, 1)')
  }
  if (length(step) != 1 || is.na(step) || step <= 0 || step >= alpha) {
    stop('step must be a single value in (0, alpha)')
  }
  psi <- sort(unique(psi))
  a <- check_paired_args(N, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative, margin)
  Test <- a$Test
  N <- a$N
  # The Berger-Boos parameter of the final analysis matters only for the unconditional test
  bb.test <- if (Test == 'McNemar-uncond') bb.gamma else 0
  n1 <- paired_n1(N, omega, n.interim)
  # The null hypothesis has one boundary, or two for a test of equivalence. A point of the
  # null hypothesis is the pair c(psi, k) of the probability psi of a discordant pair and
  # the number k of the boundary, so that the probe passed by bisect_level() carries both
  sides <- null_sides(alternative, margin, 'RD')
  grids <- lapply(sides, function(s) {
    bd <- paired_boundary(psi, s$alternative, s$margin)
    if (!any(bd$ok)) stop('no value of psi lies on the boundary of the null hypothesis')
    list(psi = psi[bd$ok], interval = bd$range)
  })
  boundary_at <- function(x, k) {
    s <- sides[[k]]
    paired_boundary(x, s$alternative, s$margin)
  }
  # Rates on the two boundaries that differ by at most tie.tol count as equal and the
  # lower boundary is kept, as in BinaryAlphaAdjBSSR()
  tie.tol <- if (length(sides) == 2) 1e-12 else 0
  better <- function(y, k, m) {
    is.null(m) || y > m$y + tie.tol || isTRUE(y >= m$y - tie.tol && k < m$x[2])
  }
  # Largest rate of a design over the grid, refined between the grid points unless
  # maximize = 'grid'. With a probe c(psi, k) the rate is first evaluated there, and a
  # value above alpha is returned at once
  max_tie <- function(des, probe = NULL) {
    if (!is.null(probe)) {
      v <- des$f(probe[1], probe[2])
      if (isTRUE(v > alpha)) return(list(x = probe, y = v, bound = NA_real_))
    }
    m <- NULL
    for (k in seq_along(grids)) {
      ps <- grids[[k]]$psi
      y <- des$f(ps, k)
      mk <- if (maximize == 'grid') {
        list(x = ps[which.max(y)], y = max(y))
      } else {
        refine_max(function(x) des$f(x, k), ps, y)
      }
      if (better(mk$y, k, m)) m <- list(x = c(mk$x, k), y = mk$y)
    }
    m$bound <- NA_real_
    m
  }
  # Certified maximum over the intervals, combined with the result m of max_tie
  certify <- if (maximize == 'certified') {
    function(des, m) {
      bound <- m$y
      for (k in seq_along(sides)) {
        s <- sides[[k]]
        cm <- paired_tie_certify(des$rr.list, des$Ns, des$W, grids[[k]]$interval,
                                 s$alternative, s$margin)
        if (better(cm$y, k, m)) {
          m$x <- c(cm$x, k)
          m$y <- cm$y
        }
        bound <- max(bound, cm$bound)
      }
      m$bound <- max(bound, m$y)
      m
    }
  }
  passes <- function(m) m$y <= alpha && (is.na(m$bound) || m$bound <= alpha)
  region <- function(pv, level) {
    rr <- pv %<<% level
    rr[is.na(rr)] <- FALSE
    rr
  }
  # Fixed design
  pv.fixed <- get_paired_pvalue(N, Test, alternative, a$n.grid, bb.test, a$ref.pvalue,
                                a$margin)
  make_trad <- function(level) {
    rr <- region(pv.fixed, level)
    f <- function(x, k) {
      b <- boundary_at(x, k)
      paired_power(rr, b$p10, b$p01)
    }
    list(f = f, rr.list = list(rr), Ns = N, W = list(rep(1, N + 1)))
  }
  res.trad <- bisect_level(make_trad, max_tie, certify, alpha, tol)
  if (res.trad$level < alpha) res.trad$level <- report_level(list(pv.fixed), res.trad$level)
  # Design with re-estimation
  design_at <- function(level.ss) {
    map <- paired_reestimate(n1, 0:n1, Delta.A, N, tar.power, restricted, alternative,
                             n.grid, bb.gamma, ss.method, ss.Test, level.ss, N.min, N.max,
                             ref.pvalue, margin, search, search.limit)
    Ns <- sort(unique(map$N))
    list(Ns = Ns, rr.id = as.integer(match(map$N, Ns) - 1L),
         pv.list = lapply(Ns, function(n) {
           get_paired_pvalue(n, Test, alternative, a$n.grid, bb.test, a$ref.pvalue,
                             a$margin)
         }),
         W = if (maximize == 'certified') paired_tie_weights(map$N, n1, Ns) else NULL)
  }
  make_bssr <- function(d, level) {
    rr.list <- lapply(d$pv.list, region, level = level)
    f <- function(x, k) {
      b <- boundary_at(x, k)
      paired_bssr_power(rr.list, d$rr.id, n1, b$p10, b$p01)
    }
    list(f = f, rr.list = rr.list, Ns = d$Ns, W = d$W)
  }
  if (adjust == 'test') {
    d <- design_at(ss.alpha)
    res.bssr <- bisect_level(function(level) make_bssr(d, level), max_tie, certify, alpha,
                             tol)
    if (res.bssr$level < alpha) res.bssr$level <- report_level(d$pv.list, res.bssr$level)
  } else {
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
  interval <- if (length(sides) == 1) {
    grids[[1]]$interval
  } else {
    iv <- t(vapply(grids, function(g) g$interval, numeric(2)))
    dimnames(iv) <- list(c('lower', 'upper'), c('from', 'to'))
    iv
  }
  out <- data.frame(
    Design = c('BSSR', 'TRAD'),
    alpha = alpha,
    max.TIE = c(res.bssr$m0$y, res.trad$m0$y),
    alpha.adj = c(res.bssr$level, res.trad$level),
    max.TIE.adj = c(res.bssr$m$y, res.trad$m$y),
    psi.adj = c(res.bssr$m$x[1], res.trad$m$x[1]),
    max.TIE.bound = c(res.bssr$m0$bound, res.trad$m0$bound),
    max.TIE.adj.bound = c(res.bssr$m$bound, res.trad$m$bound),
    stringsAsFactors = FALSE
  )
  if (length(sides) == 2) {
    out$boundary.adj <- c('lower', 'upper')[c(res.bssr$m$x[2], res.trad$m$x[2])]
  }
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'adjust') <- adjust
  attr(out, 'maximize') <- maximize
  if (maximize == 'certified') attr(out, 'interval') <- interval
  attr(out, 'nuisance') <- 'psi'
  attr(out, 'N') <- N
  attr(out, 'n.interim') <- n1
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'ss.method') <- ss.method
  attr(out, 'search') <- search
  attr(out, 'margin') <- margin
  attr(out, 'margin.scale') <- 'RD'
  attr(out, 'ref.pvalue') <- ref.pvalue
  class(out) <- c('bbssr_alphaadj', 'data.frame')
  out
}
