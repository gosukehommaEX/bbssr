#' Conditional Rejection Probabilities of a Blinded Sample Size Re-estimation Design
#'
#' Computes, under the null hypothesis of equal response probabilities, the probability
#' that a two-arm trial with a binary endpoint and blinded sample size re-estimation
#' (BSSR) rejects the null hypothesis given the pooled number of responders at the interim
#' analysis and in the second stage. These conditional rejection probabilities do not
#' depend on the common response probability, and the type I error rate is their average
#' over the binomial distribution of the two pooled counts. The type I error rate can also
#' be decomposed by the pooled number of interim responders, which shows the interim
#' outcomes that raise it.
#'
#' @inheritParams BinaryPowerBSSR
#' @param theta Optional vector of common response probabilities at which the type I
#'   error rate is decomposed by the pooled number of interim responders. Default is
#'   \code{NULL}, which skips the decomposition
#' @param margin Must be 0, the default, since only tests of superiority are covered (see
#'   Details). The argument keeps the arguments the same as those of
#'   \code{\link{BinaryTypeIErrorBSSR}}
#' @param ref.pvalue Logical. If \code{TRUE}, the maximization over the nuisance parameter
#'   of the unconditional tests is refined between the grid points, see
#'   \code{\link{BinaryRR}}. Default is \code{FALSE}. It applies to the final analysis and
#'   to the exact re-estimation
#'
#' @return An object of class \code{bbssr_crp}, a data frame with one row for every pair
#'   of pooled responder counts that the design can produce, ordered by \code{s} and then
#'   by \code{s2}, containing:
#' \describe{
#'   \item{s}{Pooled number of responders at the interim analysis}
#'   \item{s2}{Pooled number of responders in the second stage}
#'   \item{N1}{Final sample size of group 1 reached from \code{s}}
#'   \item{N2}{Final sample size of group 2 reached from \code{s}}
#'   \item{CRP}{Rejection probability given \code{s} and \code{s2}}
#'   \item{CRP.total}{Rejection probability given only the total number of responders
#'     \code{s + s2}, when the responder count of group 1 is treated as one
#'     hypergeometric count among all \code{N1 + N2} patients}
#' }
#' The attribute \code{max} is a data frame with the largest value of each of the two
#' probabilities and the first outcome, in the order of the rows, at which it is attained
#' up to rounding error, the attribute \code{exceed} gives the
#' number of outcomes at which each of them exceeds \code{alpha}, and the attribute
#' \code{reestimation} holds the final sample size for every pooled number of interim
#' responders. When \code{theta} is supplied, the attribute \code{by.s} is a data frame
#' with one row for each value of \code{theta} and of \code{s}, containing the probability
#' \code{prob.s} of \code{s}, the type I error rate \code{TIE.s} given \code{s} and the
#' contribution \code{prob.s * TIE.s} of \code{s} to the type I error rate, and the
#' attribute \code{TIE} gives the type I error rate at each value of \code{theta}.
#'
#' @details
#' Let \code{n11} and \code{n12} be the interim sample sizes of the two groups, and let
#' \code{n21(s)} and \code{n22(s)} be the second-stage sizes reached from the pooled
#' number \code{s} of interim responders. Under the null hypothesis both groups share the
#' response probability \code{theta}. Given \code{s}, the number of interim responders in
#' group 1 is hypergeometric, being the number of group 1 patients among \code{s}
#' responders drawn from \code{n11 + n12} patients. Given the pooled number \code{s2} of
#' second-stage responders, the number of second-stage responders in group 1 is an
#' independent hypergeometric count of the same kind. \code{CRP} sums the rejection region
#' of the final sample size over these two distributions and does not depend on
#' \code{theta}. The type I error rate is
#' \code{sum_s sum_s2 b(s; n11 + n12, theta) b(s2; n21(s) + n22(s), theta) CRP(s, s2)},
#' where \code{b} is the binomial probability, so it is at most \code{alpha} for every
#' \code{theta} whenever \code{CRP} is at most \code{alpha} for every outcome.
#'
#' \code{CRP.total} is the conditional size of the final test given the total number of
#' responders, with the responder count of group 1 hypergeometric among all final
#' patients. Fisher's exact test keeps it at or below \code{alpha} by construction. The
#' allocation is fixed within each stage, so the responder count of group 1 given
#' \code{s} and \code{s2} is the sum of two hypergeometric counts rather than one, and
#' \code{CRP} differs from \code{CRP.total}. When the final sample size does not depend
#' on \code{s}, averaging \code{CRP} over the hypergeometric distribution of \code{s}
#' given the total gives \code{CRP.total}. Re-estimation makes the final sample size
#' depend on \code{s}, and this averaging no longer holds. Outcomes with \code{CRP} above
#' \code{alpha} therefore occur in a fixed-sample design as well, and whether they raise
#' the type I error rate is seen from the decomposition by \code{s} rather than from
#' \code{CRP} alone.
#'
#' Only tests of superiority are covered. With a non-inferiority margin the two response
#' probabilities differ on the boundary of the null hypothesis, and the conditional
#' distributions of the responder counts then depend on them.
#'
#' @examples
#' crp <- BinaryCondRejectBSSR(
#'   Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Fisher', ss.method = 'standard',
#'   theta = c(0.3, 0.5)
#' )
#' print(crp)
#' attr(crp, 'TIE')
#'
#' \donttest{
#' crp <- BinaryCondRejectBSSR(
#'   Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Boschloo', ss.method = 'standard'
#' )
#' attr(crp, 'max')
#' plot(crp)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryTypeIErrorBSSR}}, \code{\link{BinaryPowerBSSR}}
#' @export
#' @import fpCompare
#' @importFrom stats dbinom
BinaryCondRejectBSSR <- function(Delta.A, N1, N2, omega = NULL, r, alpha, tar.power, Test,
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
                                 theta = NULL, margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  rounding <- match.arg(rounding)
  if (length(margin) != 1 || !is.numeric(margin) || is.na(margin) || margin != 0) {
    stop('only tests of superiority are covered, so margin must be 0')
  }
  if (!is.null(theta)) {
    if (length(theta) < 1 || anyNA(theta) || any(theta < 0 | theta > 1)) {
      stop('theta must be NULL or a vector of values in [0, 1]')
    }
    theta <- sort(unique(theta))
  }
  map <- bssr_map(Delta.A, N1, N2, omega, n.interim, r, alpha, tar.power, Test,
                  restricted, alternative, tsmethod, n.grid, bb.gamma, effect,
                  ss.method, ss.Test, ss.alpha, rounding, N.min, N.max, ref.pvalue, 0)
  setup <- bssr_setup(map)
  rr.list <- lapply(seq_along(setup$N1), function(k) {
    get_rr(setup$N1[k], setup$N2[k], alpha, Test, alternative, tsmethod, n.grid, bb.gamma,
           ref.pvalue, 0)
  })
  # Index of the final sample sizes reached from each pooled number of interim responders
  id.s <- match(paste(map$N1, map$N2, sep = '_'), paste(setup$N1, setup$N2, sep = '_'))
  cr <- bssr_cond_reject(rr.list, as.integer(id.s - 1L), setup$n21, setup$n22,
                         setup$n11, setup$n12)
  out <- data.frame(s = cr$s, s2 = cr$s2, N1 = map$N1[cr$s + 1L], N2 = map$N2[cr$s + 1L],
                    CRP = cr$crp, CRP.total = cr$crp.total)
  # Largest value of each probability, at the first outcome where it is attained. Values
  # that are equal in exact arithmetic can differ in the last bits, so the maximum is
  # located with a relative tolerance
  top <- function(v) {
    i <- which(v >= max(v) * (1 - 1e-12))[1]
    data.frame(value = v[i], s = out$s[i], s2 = out$s2[i], N1 = out$N1[i],
               N2 = out$N2[i])
  }
  attr(out, 'max') <- cbind(Quantity = c('CRP', 'CRP.total'),
                            rbind(top(out$CRP), top(out$CRP.total)),
                            stringsAsFactors = FALSE)
  attr(out, 'exceed') <- c(CRP = sum(out$CRP %>>% alpha),
                           CRP.total = sum(out$CRP.total %>>% alpha))
  if (!is.null(theta)) {
    n1 <- setup$n11 + setup$n12
    n2.row <- out$N1 + out$N2 - n1
    parts <- lapply(theta, function(t) {
      # Type I error rate given s, from the binomial distribution of s2
      cond <- rowsum(dbinom(out$s2, n2.row, t) * out$CRP, out$s, reorder = TRUE)[, 1]
      prob <- dbinom(map$s, n1, t)
      data.frame(theta = t, s = map$s, N1 = map$N1, N2 = map$N2, prob.s = prob,
                 TIE.s = unname(cond), contribution = prob * unname(cond))
    })
    by.s <- do.call(rbind, parts)
    rownames(by.s) <- NULL
    attr(out, 'by.s') <- by.s
    attr(out, 'TIE') <- data.frame(
      theta = theta,
      TIE = vapply(parts, function(d) sum(d$contribution), numeric(1))
    )
  }
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
  attr(out, 'ref.pvalue') <- ref.pvalue
  attr(out, 'reestimation') <- map
  class(out) <- c('bbssr_crp', 'data.frame')
  out
}
