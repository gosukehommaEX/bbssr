#' Type I Error Rate of a Blinded Sample Size Re-estimation Design
#'
#' Evaluates the type I error rate of a two-arm trial with a binary endpoint and blinded
#' sample size re-estimation (BSSR) over the common response probability, or over the
#' pooled response probability on the boundary of a non-inferiority hypothesis, together
#' with that of the corresponding fixed-sample design, and locates the largest value of
#' each.
#'
#' @inheritParams BinaryPowerBSSR
#' @param theta Grid of common response probabilities at which the type I error rate is
#'   evaluated. Default is \code{seq(0, 1, by = 0.005)}. With a non-zero
#'   \code{margin}, or a margin on the scale of the risk ratio, \code{theta} is the pooled
#'   response probability
#'   \code{(r p1 + p2) / (1 + r)} on the boundary of the null hypothesis, see Details
#' @param maximize How the largest type I error rate is located. \code{'certified'}
#'   (default) finds it over the whole interval from the smallest to the largest value of
#'   \code{theta}, restricted with a margin to the values at which both response
#'   probabilities on the null boundary lie in the unit interval, together with an upper
#'   bound, see Details. \code{'refined'} refines the largest local maxima on the grid
#'   by a one-dimensional optimization between the neighbouring grid points, and
#'   \code{'grid'} takes the largest value on the grid
#'
#' @return An object of class \code{bbssr_tie}, a data frame with one row per element of
#'   \code{theta} containing:
#' \describe{
#'   \item{theta}{Common response probability of the two groups, or the pooled response
#'     probability on the null boundary}
#'   \item{p1}{Response probability of group 1}
#'   \item{p2}{Response probability of group 2}
#'   \item{TIE.BSSR}{Type I error rate of the BSSR design}
#'   \item{TIE.TRAD}{Type I error rate of the fixed-sample design with sample sizes
#'     \code{N1} and \code{N2}}
#' }
#' The attribute \code{max} is a data frame with one row for each design, \code{'BSSR'}
#' and \code{'TRAD'} in the column \code{Design}, holding the largest type I error rate
#' \code{TIE}, the common response probability, or pooled response probability on the
#' null boundary, \code{theta} at which it occurs and, with
#' \code{maximize = 'certified'}, the upper bound \code{bound} of the type I error rate
#' over the interval given by the attribute \code{interval} (otherwise \code{NA}).
#' The attribute \code{reestimation} holds the final sample size for every pooled number
#' of interim responders.
#'
#' For \code{alternative = 'equivalence'} the null hypothesis has two boundaries, the
#' lower one \code{p1 - p2 = -lower} (\code{p1 = lower p2} on the scale of the risk
#' ratio) and the upper one \code{p1 - p2 = upper} (\code{p1 = upper p2}). The data
#' frame then has a first column \code{boundary}, \code{'lower'} or \code{'upper'},
#' and a row for each element of \code{theta} on each boundary. The rejection
#' probability is that of the two one-sided tests together. The largest rate of each
#' design in the attribute \code{max} is the larger of its maxima on the two boundaries,
#' its bound the larger of the two bounds, and the column \code{boundary} of the
#' attribute gives the boundary on which it occurs. The attribute \code{interval} is
#' then a matrix with a row for each boundary.
#'
#' @details
#' Under the null hypothesis both groups share the response probability \code{theta}. The
#' sample size is re-estimated exactly as in \code{\link{BinaryPowerBSSR}}, from the
#' assumed effect \code{Delta.A}, and the rejection probability is summed over every
#' interim outcome and every outcome of the second stage.
#'
#' The type I error rate is a polynomial in \code{theta}, since both response
#' probabilities are linear in \code{theta}. With \code{maximize = 'certified'} the
#' polynomial is expressed in the Bernstein basis over the interval from the smallest to
#' the largest value of \code{theta}, restricted with a margin to the values at which both
#' response probabilities lie in the unit interval. The basis polynomials are
#' non-negative and sum to one, so on any subinterval the polynomial does not exceed its
#' largest Bernstein coefficient there, and its first and last coefficients are its values
#' at the two ends. The interval is halved repeatedly by the algorithm of de Casteljau, and
#' a subinterval is set aside once its largest coefficient exceeds the largest value found
#' by at most 1e-12. The largest coefficient of the subintervals set aside, reported as
#' \code{bound}, is an upper bound of the type I error rate at every value of \code{theta}
#' in the interval, up to rounding error, and \code{TIE} is within 1e-12 of it. If 1e5
#' subdivisions do not suffice, a warning is given; \code{bound} remains an upper bound,
#' but \code{TIE} can then be more than 1e-12 below it. With \code{maximize = 'refined'}
#' the three largest local maxima on the grid are refined, which usually finds the maximum
#' when the grid is fine compared with the spacing of the local maxima but does not
#' guarantee it, and with \code{maximize = 'grid'} the largest value on the grid is
#' reported.
#'
#' With a non-inferiority \code{margin} the null hypothesis is \code{p1 - p2 <= -margin},
#' or \code{p1 - p2 >= margin} for \code{alternative = 'less'}, and the type I error rate
#' is evaluated on its boundary, as in Friede et al. (2007). The boundary is parametrized
#' by the pooled response probability \code{theta}, and the values of \code{theta} at
#' which a response probability falls outside the unit interval are dropped. With
#' \code{margin.scale = 'RR'} the null hypothesis is \code{p1 / p2 <= margin}, or
#' \code{p1 / p2 >= margin} for \code{alternative = 'less'}, and its boundary
#' \code{p1 = margin p2} is parametrized in the same way, with
#' \code{p2 = (1 + r) theta / (1 + r margin)}. Both response probabilities are again
#' linear in \code{theta}, so the certified maximum applies on either scale.
#'
#' @references
#' Friede T, Kieser M (2004). Sample size recalculation for binary data in internal pilot
#' study designs. \emph{Pharmaceutical Statistics}, 3(4), 269-279.
#'
#' Kieser M (2020). \emph{Methods and Applications of Sample Size Calculation and
#' Recalculation in Clinical Trials}. Springer, Cham.
#'
#' Friede T, Mitchell C, Mueller-Velten G (2007). Blinded sample size reestimation in
#' non-inferiority trials with binary endpoints. \emph{Biometrical Journal}, 49(6),
#' 903-916.
#'
#' @examples
#' tie <- BinaryTypeIErrorBSSR(
#'   Delta.A = 0.3, N1 = 20, N2 = 20, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
#'   theta = seq(0.05, 0.95, by = 0.05)
#' )
#' print(tie)
#'
#' \donttest{
#' tie <- BinaryTypeIErrorBSSR(
#'   Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Boschloo'
#' )
#' attr(tie, 'max')
#' plot(tie)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryAlphaAdjBSSR}}, \code{\link{BinaryPowerBSSR}}
#' @export
#' @importFrom stats dbinom
BinaryTypeIErrorBSSR <- function(Delta.A, N1, N2, omega = NULL, r, alpha, tar.power, Test,
                                 restricted = FALSE,
                                 alternative = c('greater', 'less', 'two.sided', 'equivalence'),
                                 tsmethod = c('minlike', 'central', 'blaker'),
                                 n.grid = 100, bb.gamma = 0,
                                 effect = c('RD', 'RR', 'OR'),
                                 ss.method = c('exact', 'standard', 'null.variance',
                                               'alternative.variance'),
                                 ss.Test = Test, ss.alpha = alpha,
                                 rounding = c('group', 'friede-kieser', 'total',
                                              'nearest'),
                                 search = c('crossing', 'smallest', 'stable'),
                                 search.limit = c(2, 50),
                                 N.min = NULL, N.max = NULL, n.interim = NULL,
                                 theta = seq(0, 1, by = 0.005),
                                 maximize = c('certified', 'refined', 'grid'),
                                 margin = 0, ref.pvalue = FALSE,
                                 margin.scale = c('RD', 'RR')) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  search <- match.arg(search)
  maximize <- match.arg(maximize)
  margin.scale <- match.arg(margin.scale)
  margin <- norm_margin(alternative, margin, margin.scale)
  if (length(theta) < 1 || anyNA(theta) || any(theta < 0 | theta > 1)) {
    stop('theta must be a vector of values in [0, 1]')
  }
  theta <- sort(unique(theta))
  map <- bssr_map(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                  restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                  ss.method, ss.Test, ss.alpha, rounding, N.min, N.max, ref.pvalue,
                  margin, search = search, search.limit = search.limit,
                  margin.scale = margin.scale)
  # The null hypothesis has one boundary, or two for a test of equivalence. On each the
  # response probabilities are functions of the pooled probability theta, both equal to
  # theta when the margin is 0
  sides <- null_sides(alternative, margin, margin.scale)
  for (s in sides) {
    if (!any(null_boundary(theta, r, s$alternative, s$margin, margin.scale)$ok)) {
      stop('no value of theta gives response probabilities in [0, 1] on the null boundary')
    }
  }
  setup <- bssr_setup(map)
  rr.list <- lapply(seq_along(setup$N1), function(k) {
    get_rr(setup$N1[k], setup$N2[k], alpha, Test, alternative, tsmethod, n.grid, bb.gamma,
           ref.pvalue, margin, margin.scale)
  })
  rr.fixed <- get_rr(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma,
                     ref.pvalue, margin, margin.scale)
  # Weights of the certified maximization, which do not depend on the boundary. The
  # fixed-sample design is certified as a design without an interim stage
  w.bssr <- w.trad <- NULL
  if (maximize == 'certified') {
    w.bssr <- tie_weights(setup, rr.list)
    w.trad <- tie_weights(fixed_setup(N1, N2), list(rr.fixed))
  }
  per.side <- lapply(sides, function(s) {
    # Interval over which the largest rate is certified, and the grid on the boundary
    interval <- null_range(theta, r, s$alternative, s$margin, margin.scale)
    th <- theta[null_boundary(theta, r, s$alternative, s$margin, margin.scale)$ok]
    f.bssr <- function(t) {
      b <- null_boundary(t, r, s$alternative, s$margin, margin.scale)
      bssr_reject(setup, rr.list, b$p1, b$p2)
    }
    f.trad <- function(t) {
      b <- null_boundary(t, r, s$alternative, s$margin, margin.scale)
      vapply(seq_along(t), function(i) {
        power_from_rr(rr.fixed, dbinom(0:N1, N1, b$p1[i]), dbinom(0:N2, N2, b$p2[i]))
      }, numeric(1))
    }
    tie.bssr <- f.bssr(th)
    tie.trad <- f.trad(th)
    # Largest rate on the grid, refined between the grid points or certified over the
    # interval
    locate <- function(f, y, w) {
      m <- if (maximize == 'refined') {
        refine_max(f, th, y)
      } else {
        list(x = th[which.max(y)], y = max(y))
      }
      m$bound <- NA_real_
      if (maximize == 'certified') {
        cm <- tie_certify(w, interval, r, s$alternative, s$margin,
                          margin.scale = margin.scale)
        if (cm$y > m$y) {
          m$x <- cm$x
          m$y <- cm$y
        }
        m$bound <- max(cm$bound, m$y)
      }
      m
    }
    b <- null_boundary(th, r, s$alternative, s$margin, margin.scale)
    list(df = data.frame(theta = th, p1 = b$p1, p2 = b$p2, TIE.BSSR = tie.bssr,
                         TIE.TRAD = tie.trad),
         m.bssr = locate(f.bssr, tie.bssr, w.bssr),
         m.trad = locate(f.trad, tie.trad, w.trad),
         interval = interval, side = s$side)
  })
  # The largest rate of a design is the larger of its maxima on the boundaries, and the
  # bound the larger of their bounds. Maxima that differ by at most 1e-12 count as equal
  # and the lower boundary is reported, so that the result for a symmetric design does
  # not depend on rounding
  pick <- function(d) {
    y <- vapply(per.side, function(p) p[[d]]$y, numeric(1))
    k <- if (length(y) == 1) 1L else which(y >= max(y) - 1e-12)[1]
    m <- per.side[[k]][[d]]
    m$bound <- max(vapply(per.side, function(p) p[[d]]$bound, numeric(1)))
    m$side <- per.side[[k]]$side
    m
  }
  m.bssr <- pick('m.bssr')
  m.trad <- pick('m.trad')
  mx <- data.frame(Design = c('BSSR', 'TRAD'),
                   theta = c(m.bssr$x, m.trad$x),
                   TIE = c(m.bssr$y, m.trad$y),
                   bound = c(m.bssr$bound, m.trad$bound),
                   stringsAsFactors = FALSE)
  if (length(per.side) == 1) {
    out <- per.side[[1]]$df
    interval <- per.side[[1]]$interval
  } else {
    # For a test of equivalence the rates on the two boundaries are stacked
    out <- do.call(rbind, lapply(per.side, function(p) {
      data.frame(boundary = p$side, p$df, stringsAsFactors = FALSE)
    }))
    rownames(out) <- NULL
    interval <- t(vapply(per.side, function(p) p$interval, numeric(2)))
    dimnames(interval) <- list(c('lower', 'upper'), c('from', 'to'))
    mx$boundary <- c(m.bssr$side, m.trad$side)
  }
  attr(out, 'max') <- mx
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N1') <- N1
  attr(out, 'N2') <- N2
  attr(out, 'n1.interim') <- setup$n11
  attr(out, 'n2.interim') <- setup$n12
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'search') <- search
  attr(out, 'maximize') <- maximize
  if (maximize == 'certified') attr(out, 'interval') <- interval
  attr(out, 'margin') <- margin
  attr(out, 'margin.scale') <- margin.scale
  attr(out, 'ref.pvalue') <- ref.pvalue
  attr(out, 'reestimation') <- map
  class(out) <- c('bbssr_tie', 'data.frame')
  out
}
