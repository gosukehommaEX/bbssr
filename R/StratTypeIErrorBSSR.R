#' Type I Error Rate of Blinded Sample Size Re-Estimation for a Stratified Design
#'
#' Calculates the type I error rate of a two-arm trial with a binary endpoint, stratified
#' into two or more strata, with blinded sample size re-estimation, over the common
#' response probabilities of the strata under the null hypothesis, together with that of
#' the fixed design with the planned group sizes. With two strata the rate is computed
#' exactly and its largest value over a rectangle of probabilities is certified; with three
#' or more strata the rate of the re-estimation design is estimated by simulation.
#'
#' @param theta Common response probabilities of the strata under the null hypothesis. A
#'   vector gives, with two strata, the grid of all pairs of its values, and with three or
#'   more strata the points at which every stratum has the same value. A matrix with one
#'   column per stratum gives the points directly. Default is
#'   \code{seq(0, 1, by = 0.05)}
#' @param maximize How the largest type I error rate is located. With
#'   \code{'certified'} (the default with two strata, and available only with two) it is
#'   certified over the rectangle from the smallest to the largest value of \code{theta} in
#'   each stratum, see Details. With \code{'grid'} (the default with three or more strata)
#'   it is the largest value at the points of \code{theta}
#' @param nsim Number of simulated trials at each point with three or more strata. Default
#'   is \code{1e5}. Ignored with two strata
#' @inheritParams StratPowerBSSR
#'
#' @return An object of class \code{bbssr_stie}, a data frame with one row per point,
#'   containing the probabilities \code{theta.S1, theta.S2, ...}, the type I error rates
#'   \code{TIE.BSSR} and \code{TIE.TRAD} of the designs with re-estimation and without it,
#'   and, with three or more strata or a random allocation, the Monte Carlo standard error
#'   \code{se.BSSR}, and with three or more strata and a random allocation \code{se.TRAD}
#'   for \code{TIE.TRAD}. The attribute \code{max} is a data frame with the location, the
#'   value and, with \code{maximize = 'certified'}, the upper bound \code{bound} of the
#'   largest rate of each design, and with a random allocation the standard error
#'   \code{se} of the rate at the point of \code{theta} with the largest rate (\code{NA}
#'   for the exact rate of the fixed design with two strata).
#'
#' @details
#' Under the null hypothesis both groups of stratum \code{j} respond with the probability
#' \code{theta_j}. Given the interim and the second-stage numbers of responders of a
#' stratum, the numbers of responders of group 1 are hypergeometric whatever
#' \code{theta_j}, so with two strata the type I error rate is a polynomial in
#' \code{(theta_1, theta_2)}. Its coefficients in the products of the Bernstein bases of the
#' two strata are computed exactly, and the rate is evaluated at the points of
#' \code{theta} from them. Over a rectangle the polynomial lies between the smallest and the
#' largest of its coefficients there, so with \code{maximize = 'certified'} the rectangle is
#' halved along its longer side until the largest coefficient of every piece exceeds the
#' largest value found by at most 1e-12. The largest value is reported with this upper
#' bound, which the type I error rate cannot exceed anywhere in the rectangle. The fixed
#' design is certified in the same way.
#'
#' With three or more strata the type I error rate of the re-estimation design is estimated
#' at each point from \code{nsim} simulated trials, with the random numbers of R, and that
#' of the fixed design with the fixed allocation is computed exactly.
#'
#' With a random allocation and two strata the allocation sequences are drawn as in
#' \code{\link{StratPowerBSSR}}. Given a sequence the type I error rate is a polynomial
#' as above, and the rate of the re-estimation design is the average of these polynomials
#' over \code{nseq} sequences, itself a polynomial, whose largest value is certified as
#' above. The certificate holds for this average, and \code{se.BSSR} gives the Monte Carlo
#' standard error of the average at each point, from the variation between the
#' sequences; the attribute \code{max} gives it at the point of \code{theta} with the
#' largest rate. The rate of the fixed design with the planned total and the same
#' allocation is exact: the sum over the stratum and group sizes of their probability
#' times the polynomial of the fixed design they give.
#'
#' With a random allocation and three or more strata (\code{alloc = 'stratum'} or
#' \code{'both'}) the allocation is simulated with the responses as in
#' \code{\link{StratPowerBSSR}}, and the rates of both designs are estimated at each point
#' from \code{nsim} trials each.
#'
#' @examples
#' tie <- StratTypeIErrorBSSR(Delta.A = 0.2, N1 = c(20, 20), N2 = c(20, 20), omega = 0.5,
#'                            r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                            ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.2))
#' attr(tie, 'max')
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratAlphaAdjBSSR}}, \code{\link{StratPowerBSSR}}
#' @export
StratTypeIErrorBSSR <- function(Delta.A, N1, N2, omega = NULL, r, alpha, tar.power, Test,
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
                                maximize = c('certified', 'grid'), nsim = 1e5,
                                alloc = c('fixed', 'stratum', 'group', 'both'),
                                block = 4, ss.frac = c('observed', 'planned'),
                                nseq = 100) {
  # Without an explicit choice, 'grid' is used with three or more strata
  maximize.given <- !missing(maximize)
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  ss.strata <- match.arg(ss.strata)
  maximize <- match.arg(maximize)
  alloc <- match.arg(alloc)
  ss.frac <- match.arg(ss.frac)
  d <- strat_design(N1, N2, r, frac, omega, n.interim, alpha, Test, alternative, tsmethod,
                    ss.method, nsim)
  a <- d$a
  J <- d$J
  if (alloc != 'fixed' && J > 2) {
    # The allocation is simulated with the responses
    st <- list(al = check_strat_alloc(alloc, NULL, d$frac, d$r, block, J),
               N.plan = sum(a$N1 + a$N2), n = sum(d$n1 + d$n2))
    if (!is.null(N.max) && N.max < st$n) {
      stop('N.max must be at least the interim total sample size')
    }
  }
  if (J > 2 && !maximize.given) maximize <- 'grid'
  if (J > 2 && maximize == 'certified') {
    stop("maximize = 'certified' requires two strata; use maximize = 'grid'")
  }
  pts <- strat_theta(theta, J)
  box <- apply(pts, 2, range)
  if (maximize == 'certified' && any(box[2, ] - box[1, ] <= 0)) {
    stop("maximize = 'certified' requires theta to span an interval in each stratum")
  }
  re.args <- list(Delta.A = Delta.A, effect = effect, frac = d$frac, r = d$r,
                  ss.method = ss.method, ss.Test = ss.Test, ss.alpha = ss.alpha,
                  tar.power = tar.power, N1.plan = a$N1, N2.plan = a$N2,
                  restricted = restricted, N.min = N.min, N.max = N.max, search = search,
                  search.limit = search.limit, ss.strata = ss.strata)
  se.BSSR <- NULL
  se.TRAD <- NULL
  map <- NULL
  pick <- function(y) {
    k <- which.max(y)
    list(x = pts[k, ], y = y[k], bound = NA_real_)
  }
  if (J == 2) {
    if (alloc == 'fixed') {
      sm <- strat_map(d$n1, d$n2, alpha, a$Test, a$alternative, a$tsmethod, re.args)
      map <- sm$map
      coef.bssr <- strat_tie_coef(d$n1, d$n2, sm$designs, sm$id)
      crit <- get_strat_crit(a$N1, a$N2, alpha, a$Test, a$alternative, a$tsmethod)
      coef.trad <- strat_tie_coef(c(0L, 0L), c(0L, 0L), list(crit), 0L)
    } else {
      st <- strat_alloc_setup(d, alloc, block, nseq, function(level) re.args, ss.frac,
                              restricted, N.min, N.max)
      lay <- st$layout(ss.alpha)
      res <- strat_seq_tie_coef(lay$inputs, a$type, a$alt, a$ts, alpha, pts[, 1],
                                pts[, 2])
      coef.bssr <- res$coef
      se.BSSR <- apply(res$values, 2, stats::sd) / sqrt(nseq)
      coef.trad <- strat_alloc_tie_coef(st$N.plan, st$al, alloc, alpha, a$Test,
                                        a$alternative, a$tsmethod)$coef
    }
    tie.bssr <- strat_bern_eval(coef.bssr, pts[, 1], pts[, 2])
    tie.trad <- strat_bern_eval(coef.trad, pts[, 1], pts[, 2])
    locate <- function(coef, y) {
      m <- pick(y)
      if (maximize == 'certified') {
        cm <- strat_tie_certify(coef, box)
        if (cm$y > m$y) {
          m$x <- cm$x
          m$y <- cm$y
        }
        m$bound <- max(cm$bound, m$y)
      }
      m
    }
    m.bssr <- locate(coef.bssr, tie.bssr)
    m.trad <- locate(coef.trad, tie.trad)
  } else {
    sims <- lapply(seq_len(nrow(pts)), function(k) {
      if (alloc == 'fixed') {
        strat_simulate(d$n1, d$n2, pts[k, ], pts[k, ], d$nsim, alpha, a, re.args)
      } else {
        strat_alloc_simulate(st$n, st$N.plan, st$al, alloc, pts[k, ], pts[k, ], d$nsim,
                             alpha, a, re.args, ss.frac, restricted, N.min, N.max)
      }
    })
    tie.bssr <- vapply(sims, function(s) mean(s$reject), numeric(1))
    se.BSSR <- sqrt(tie.bssr * (1 - tie.bssr) / d$nsim)
    if (alloc == 'fixed') {
      tie.trad <- strat_fixed_power(a$N1, a$N2, pts, pts, alpha, a$Test, a$alternative,
                                    a$tsmethod)
    } else {
      tie.trad <- vapply(sims, function(s) mean(s$reject.trad), numeric(1))
      se.TRAD <- sqrt(tie.trad * (1 - tie.trad) / d$nsim)
    }
    m.bssr <- pick(tie.bssr)
    m.trad <- pick(tie.trad)
  }
  out <- data.frame(pts, TIE.BSSR = tie.bssr, TIE.TRAD = tie.trad)
  if (!is.null(se.BSSR)) out$se.BSSR <- se.BSSR
  if (!is.null(se.TRAD)) out$se.TRAD <- se.TRAD
  mx <- data.frame(Design = c('BSSR', 'TRAD'), rbind(m.bssr$x, m.trad$x),
                   TIE = c(m.bssr$y, m.trad$y), bound = c(m.bssr$bound, m.trad$bound),
                   stringsAsFactors = FALSE)
  names(mx)[1 + seq_len(J)] <- colnames(pts)
  if (alloc != 'fixed') {
    mx$se <- c(se.BSSR[which.max(tie.bssr)],
               if (is.null(se.TRAD)) NA_real_ else se.TRAD[which.max(tie.trad)])
  }
  attr(out, 'max') <- mx
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'tsmethod') <- a$tsmethod
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N1') <- a$N1
  attr(out, 'N2') <- a$N2
  attr(out, 'n1.interim') <- d$n1
  attr(out, 'n2.interim') <- d$n2
  attr(out, 'frac') <- d$frac
  attr(out, 'r') <- d$r
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.strata') <- ss.strata
  attr(out, 'search') <- search
  attr(out, 'maximize') <- maximize
  if (maximize == 'certified') {
    attr(out, 'interval') <- t(box)
  }
  attr(out, 'nsim') <- if (J > 2) d$nsim else NULL
  attr(out, 'reestimation') <- map
  attr(out, 'alloc') <- alloc
  if (alloc != 'fixed') {
    attr(out, 'block') <- st$al$block
    attr(out, 'ss.frac') <- ss.frac
    if (J == 2) attr(out, 'nseq') <- as.integer(nseq)
    attr(out, 'n.interim') <- as.integer(st$n)
  }
  class(out) <- c('bbssr_stie', 'data.frame')
  out
}
