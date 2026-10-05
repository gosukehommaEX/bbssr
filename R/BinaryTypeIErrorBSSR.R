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
#'   evaluated. Default is \code{seq(0.005, 0.995, by = 0.005)}. With a non-zero
#'   \code{margin}, \code{theta} is the pooled response probability
#'   \code{(r p1 + p2) / (1 + r)} on the boundary of the null hypothesis, see Details
#' @param refine Logical. If \code{TRUE} (default), the largest local maxima on the grid
#'   are refined by a one-dimensional optimization between the neighbouring grid points
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
#' The attribute \code{max} is a data frame with the largest type I error rate of each
#' design and the common response probability at which it occurs, and the attribute
#' \code{reestimation} holds the final sample size for every pooled number of interim
#' responders.
#'
#' @details
#' Under the null hypothesis both groups share the response probability \code{theta}. The
#' sample size is re-estimated exactly as in \code{\link{BinaryPowerBSSR}}, from the
#' assumed effect \code{Delta.A}, and the rejection probability is summed over every
#' interim outcome and every outcome of the second stage. Exact tests control the type I
#' error rate of a fixed-sample design, but this property is not inherited by a design in
#' which the final sample size depends on the interim data, so the rate is worth checking
#' for every design.
#'
#' The type I error rate is a polynomial in \code{theta}. With \code{refine = TRUE} the
#' three largest local maxima on the grid are refined, so the reported maximum does not
#' depend on the spacing of the grid as long as the grid separates the local maxima.
#'
#' With a non-inferiority \code{margin} the null hypothesis is \code{p1 - p2 <= -margin},
#' or \code{p1 - p2 >= margin} for \code{alternative = 'less'}, and the type I error rate
#' is evaluated on its boundary, as in Friede et al. (2007). The boundary is parametrized
#' by the pooled response probability \code{theta}, and the values of \code{theta} at
#' which a response probability falls outside the unit interval are dropped.
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
                                 alternative = c('greater', 'less', 'two.sided'),
                                 tsmethod = c('minlike', 'central'),
                                 n.grid = 100, bb.gamma = 0,
                                 effect = c('RD', 'RR', 'OR'),
                                 ss.method = c('exact', 'standard', 'null.variance',
                                               'alternative.variance'),
                                 ss.Test = Test, ss.alpha = alpha,
                                 rounding = c('group', 'friede-kieser', 'total',
                                              'nearest'),
                                 N.min = NULL, N.max = NULL, n.interim = NULL,
                                 theta = seq(0.005, 0.995, by = 0.005), refine = TRUE,
                                 margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  if (length(theta) < 1 || anyNA(theta) || any(theta < 0 | theta > 1)) {
    stop('theta must be a vector of values in [0, 1]')
  }
  theta <- sort(unique(theta))
  map <- bssr_map(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                  restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                  ss.method, ss.Test, ss.alpha, rounding, N.min, N.max, ref.pvalue,
                  margin)
  # Response probabilities on the boundary of the null hypothesis, which are both theta
  # when the margin is 0
  if (!any(null_boundary(theta, r, alternative, margin)$ok)) {
    stop('no value of theta gives response probabilities in [0, 1] on the null boundary')
  }
  theta <- theta[null_boundary(theta, r, alternative, margin)$ok]
  setup <- bssr_setup(map)
  rr.list <- lapply(seq_along(setup$N1), function(k) {
    get_rr(setup$N1[k], setup$N2[k], alpha, Test, alternative, tsmethod, n.grid, bb.gamma,
           ref.pvalue, margin)
  })
  rr.fixed <- get_rr(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma,
                     ref.pvalue, margin)
  f.bssr <- function(t) {
    b <- null_boundary(t, r, alternative, margin)
    bssr_reject(setup, rr.list, b$p1, b$p2)
  }
  f.trad <- function(t) {
    b <- null_boundary(t, r, alternative, margin)
    vapply(seq_along(t), function(i) {
      power_from_rr(rr.fixed, dbinom(0:N1, N1, b$p1[i]), dbinom(0:N2, N2, b$p2[i]))
    }, numeric(1))
  }
  tie.bssr <- f.bssr(theta)
  tie.trad <- f.trad(theta)
  if (refine) {
    m.bssr <- refine_max(f.bssr, theta, tie.bssr)
    m.trad <- refine_max(f.trad, theta, tie.trad)
  } else {
    m.bssr <- list(x = theta[which.max(tie.bssr)], y = max(tie.bssr))
    m.trad <- list(x = theta[which.max(tie.trad)], y = max(tie.trad))
  }
  b <- null_boundary(theta, r, alternative, margin)
  out <- data.frame(theta = theta, p1 = b$p1, p2 = b$p2, TIE.BSSR = tie.bssr,
                    TIE.TRAD = tie.trad)
  attr(out, 'max') <- data.frame(Design = c('BSSR', 'TRAD'),
                                 theta = c(m.bssr$x, m.trad$x),
                                 TIE = c(m.bssr$y, m.trad$y),
                                 stringsAsFactors = FALSE)
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
  attr(out, 'refine') <- refine
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  attr(out, 'reestimation') <- map
  class(out) <- c('bbssr_tie', 'data.frame')
  out
}
