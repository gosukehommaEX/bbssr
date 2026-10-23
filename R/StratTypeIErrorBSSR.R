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
#'   \code{'certified'} (default, two strata only) it is certified over the rectangle from
#'   the smallest to the largest value of \code{theta} in each stratum, see Details. With
#'   \code{'grid'} it is the largest value at the points of \code{theta}
#' @param nsim Number of simulated trials at each point with three or more strata. Default
#'   is \code{1e5}. Ignored with two strata
#' @inheritParams StratPowerBSSR
#'
#' @return An object of class \code{bbssr_stie}, a data frame with one row per point,
#'   containing the probabilities \code{theta.S1, theta.S2, ...}, the type I error rates
#'   \code{TIE.BSSR} and \code{TIE.TRAD} of the designs with re-estimation and without it,
#'   and, with three or more strata, the Monte Carlo standard error \code{se.BSSR}. The
#'   attribute \code{max} is a data frame with the location, the value and, with
#'   \code{maximize = 'certified'}, the upper bound \code{bound} of the largest rate of
#'   each design.
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
#' of the fixed design is computed exactly.
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
                                search = c('crossing', 'smallest', 'stable'),
                                search.limit = c(2, 50), N.min = NULL, N.max = NULL,
                                n.interim = NULL, frac = NULL,
                                theta = seq(0, 1, by = 0.05),
                                maximize = c('certified', 'grid'), nsim = 1e5) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  maximize <- match.arg(maximize)
  d <- strat_design(N1, N2, r, frac, omega, n.interim, alpha, Test, alternative, tsmethod,
                    ss.method, nsim)
  a <- d$a
  J <- d$J
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
                  search.limit = search.limit)
  se.BSSR <- NULL
  map <- NULL
  pick <- function(y) {
    k <- which.max(y)
    list(x = pts[k, ], y = y[k], bound = NA_real_)
  }
  if (J == 2) {
    sm <- strat_map(d$n1, d$n2, alpha, a$Test, a$alternative, a$tsmethod, re.args)
    map <- sm$map
    coef.bssr <- strat_tie_coef(d$n1, d$n2, sm$designs, sm$id)
    crit <- get_strat_crit(a$N1, a$N2, alpha, a$Test, a$alternative, a$tsmethod)
    coef.trad <- strat_tie_coef(c(0L, 0L), c(0L, 0L), list(crit), 0L)
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
      strat_simulate(d$n1, d$n2, pts[k, ], pts[k, ], d$nsim, alpha, a, re.args)
    })
    tie.bssr <- vapply(sims, function(s) mean(s$reject), numeric(1))
    se.BSSR <- sqrt(tie.bssr * (1 - tie.bssr) / d$nsim)
    tie.trad <- strat_fixed_power(a$N1, a$N2, pts, pts, alpha, a$Test, a$alternative,
                                  a$tsmethod)
    m.bssr <- pick(tie.bssr)
    m.trad <- pick(tie.trad)
  }
  out <- data.frame(pts, TIE.BSSR = tie.bssr, TIE.TRAD = tie.trad)
  if (!is.null(se.BSSR)) out$se.BSSR <- se.BSSR
  mx <- data.frame(Design = c('BSSR', 'TRAD'), rbind(m.bssr$x, m.trad$x),
                   TIE = c(m.bssr$y, m.trad$y), bound = c(m.bssr$bound, m.trad$bound),
                   stringsAsFactors = FALSE)
  names(mx)[1 + seq_len(J)] <- colnames(pts)
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
  attr(out, 'search') <- search
  attr(out, 'maximize') <- maximize
  if (maximize == 'certified') {
    attr(out, 'interval') <- t(box)
  }
  attr(out, 'nsim') <- if (J > 2) d$nsim else NULL
  attr(out, 'reestimation') <- map
  class(out) <- c('bbssr_stie', 'data.frame')
  out
}
