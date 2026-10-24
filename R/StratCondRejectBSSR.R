#' Conditional Rejection Probabilities of Blinded Sample Size Re-Estimation for a
#' Stratified Design
#'
#' Computes, under the null hypothesis of equal response probabilities in every stratum,
#' the probability that a trial with two strata and blinded sample size re-estimation
#' rejects the null hypothesis given the pooled numbers of responders of the strata at the
#' interim analysis and in the second stage. These conditional rejection probabilities do
#' not depend on the response probabilities of the strata, and the type I error rate is
#' their average over the binomial distributions of the pooled counts. With \code{theta},
#' the type I error rate is also returned for each interim outcome. This is the
#' counterpart of \code{\link{BinaryCondRejectBSSR}} for a stratified design.
#'
#' @inheritParams StratTypeIErrorBSSR
#' @param theta Optional common response probabilities of the strata at which the type I
#'   error rate is decomposed by the interim outcome: a vector, whose values form a grid
#'   of pairs, or a matrix with two columns that gives the points directly. Default is
#'   \code{NULL}, which skips the decomposition. The time of the decomposition grows with
#'   the number of points
#' @param full Logical. If \code{TRUE}, the attribute \code{full} lists the conditional
#'   rejection probability of every interim and second-stage outcome. Default is
#'   \code{FALSE}. A listing of more than 1e6 rows is refused
#'
#' @return An object of class \code{bbssr_scrp}, a data frame with one row for every pair
#'   of interim pooled numbers of responders of the strata, the first stratum varying
#'   fastest, containing:
#' \describe{
#'   \item{s1, s2}{Pooled numbers of responders of the two strata at the interim analysis}
#'   \item{N1.S1, N1.S2, N2.S1, N2.S2}{Final group sizes reached from \code{(s1, s2)}}
#'   \item{CRP.max, CRP.min}{Largest and smallest rejection probability over the pooled
#'     numbers of responders \code{(w1, w2)} of the strata in the second stage}
#'   \item{exceed}{Number of pairs \code{(w1, w2)} at which the rejection probability
#'     exceeds \code{alpha}}
#'   \item{CRP.total.max}{Largest conditional size of the final test given the final
#'     numbers of responders of the strata, see Details}
#' }
#' The attribute \code{max} is a data frame with the largest values of \code{CRP} and
#' \code{CRP.total} over all outcomes, and the first outcome, in the order of the rows and
#' then of \code{(w1, w2)} with \code{w1} varying fastest, at which each is attained up to
#' rounding error. The attribute \code{exceed} gives the numbers of outcomes
#' \code{(s1, s2, w1, w2)} at which each of the two exceeds \code{alpha}, and the
#' attribute \code{reestimation} holds the final group sizes for every interim outcome, as
#' in \code{\link{StratPowerBSSR}}. When \code{theta} is supplied, the attribute
#' \code{by.s} is a data frame with one row for each point and interim outcome,
#' containing the probability \code{prob.s} of \code{(s1, s2)}, the type I error rate
#' \code{TIE.s} given \code{(s1, s2)} and their product \code{contribution}, and the
#' attribute \code{TIE} gives the type I error rate at each point. With
#' \code{full = TRUE} the attribute \code{full} is a data frame with one row for every
#' outcome \code{(s1, s2, w1, w2)}, ordered by \code{s2}, \code{s1}, \code{w2} and
#' \code{w1}, with the columns \code{CRP} and \code{CRP.total}.
#'
#' @details
#' Let \code{n1[j]} and \code{n2[j]} be the interim group sizes of stratum \code{j}, and
#' \code{m1[j]} and \code{m2[j]} the second-stage group sizes reached from the interim
#' outcome \code{s = (s1, s2)}. Under the null hypothesis both groups of stratum \code{j}
#' respond with the probability \code{theta_j}. Given \code{s_j}, the number of interim
#' responders of group 1 in stratum \code{j} is hypergeometric, being the number of group
#' 1 patients among \code{s_j} responders drawn from the \code{n1[j] + n2[j]} interim
#' patients of the stratum. Given the pooled number \code{w_j} of second-stage responders
#' of the stratum, the number of second-stage responders of group 1 is an independent
#' hypergeometric count of the same kind. \code{CRP} sums the rejection region of the final
#' design over these four distributions and does not depend on \code{theta}. The type I
#' error rate is the sum over \code{s} and \code{w = (w1, w2)} of
#' \code{b(s1; n1[1] + n2[1], theta_1) b(s2; n1[2] + n2[2], theta_2)
#' b(w1; m1[1] + m2[1], theta_1) b(w2; m1[2] + m2[2], theta_2) CRP(s, w)}, where \code{b}
#' is the binomial probability.
#'
#' \code{CRP.total} is the conditional size of the final test given the final numbers of
#' responders \code{z_j = s_j + w_j} of the strata, with the number of responders of group
#' 1 of each stratum hypergeometric among all its final patients. The sum of the
#' contributions at a point equals the type I error rate of
#' \code{\link{StratTypeIErrorBSSR}} there.
#'
#' Only designs with two strata are covered, as for the exact computations of
#' \code{\link{StratTypeIErrorBSSR}}.
#'
#' @examples
#' crp <- StratCondRejectBSSR(Delta.A = 0.3, N1 = c(8, 8), N2 = c(8, 8), omega = 0.5,
#'                            r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                            ss.method = 'standard', theta = c(0.3, 0.5))
#' print(crp)
#' attr(crp, 'TIE')
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratTypeIErrorBSSR}}, \code{\link{BinaryCondRejectBSSR}}
#' @export
#' @importFrom stats dbinom
StratCondRejectBSSR <- function(Delta.A, N1, N2, omega = NULL, r, alpha, tar.power, Test,
                                restricted = FALSE,
                                alternative = c('greater', 'less', 'two.sided'),
                                tsmethod = c('minlike', 'central', 'blaker'),
                                effect = c('RD', 'RR', 'OR'),
                                ss.method = c('exact', 'standard'),
                                ss.Test = Test, ss.alpha = alpha,
                                ss.strata = c('stratified', 'pooled'),
                                search = c('crossing', 'smallest', 'stable'),
                                search.limit = c(2, 50), N.min = NULL, N.max = NULL,
                                n.interim = NULL, frac = NULL, theta = NULL,
                                full = FALSE) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  ss.strata <- match.arg(ss.strata)
  search <- match.arg(search)
  if (!is.logical(full) || length(full) != 1 || is.na(full)) {
    stop('full must be TRUE or FALSE')
  }
  if (length(N1) != 2 || length(N2) != 2) {
    stop('the conditional rejection probabilities are available for designs with two ',
         'strata')
  }
  d <- strat_design(N1, N2, r, frac, omega, n.interim, alpha, Test, alternative, tsmethod,
                    ss.method, 1e5)
  a <- d$a
  pts <- if (is.null(theta)) matrix(numeric(0), 0, 2) else strat_theta(theta, 2)
  re.args <- list(Delta.A = Delta.A, effect = effect, frac = d$frac, r = d$r,
                  ss.method = ss.method, ss.Test = ss.Test, ss.alpha = ss.alpha,
                  tar.power = tar.power, N1.plan = a$N1, N2.plan = a$N2,
                  restricted = restricted, N.min = N.min, N.max = N.max, search = search,
                  search.limit = search.limit, ss.strata = ss.strata)
  sm <- strat_map(d$n1, d$n2, alpha, a$Test, a$alternative, a$tsmethod, re.args)
  map <- sm$map
  if (full) {
    # Number of second-stage outcomes of every interim outcome
    rows <- sum((map$N1.S1 + map$N2.S1 - d$n1[1] - d$n2[1] + 1) *
                  (map$N1.S2 + map$N2.S2 - d$n1[2] - d$n2[2] + 1))
    if (rows > 1e6) {
      stop('the full listing would have ', format(rows), ' rows, more than 1e6; use ',
           'full = FALSE')
    }
  }
  cr <- strat_cond_reject(d$n1, d$n2, sm$designs, sm$id, alpha, full, pts[, 1], pts[, 2])
  out <- data.frame(s1 = map$s1, s2 = map$s2, N1.S1 = map$N1.S1, N1.S2 = map$N1.S2,
                    N2.S1 = map$N2.S1, N2.S2 = map$N2.S2, CRP.max = cr$crp.max,
                    CRP.min = cr$crp.min, exceed = cr$exceed,
                    CRP.total.max = cr$total.max)
  # Largest value over all outcomes, at the first interim outcome in the order of the rows
  # that attains it up to a relative 1e-12, with the first (w1, w2) of that outcome. Values
  # that are equal in exact arithmetic can differ in the last bits
  top <- function(v, w1, w2) {
    i <- which(v >= max(v) * (1 - 1e-12))[1]
    data.frame(value = v[i], s1 = out$s1[i], s2 = out$s2[i], w1 = w1[i], w2 = w2[i])
  }
  attr(out, 'max') <- cbind(Quantity = c('CRP', 'CRP.total'),
                            rbind(top(cr$crp.max, cr$w1.max, cr$w2.max),
                                  top(cr$total.max, cr$w1.total, cr$w2.total)),
                            stringsAsFactors = FALSE)
  attr(out, 'exceed') <- c(CRP = sum(cr$exceed), CRP.total = sum(cr$exceed.total))
  if (nrow(pts) > 0) {
    n.s <- d$n1 + d$n2
    parts <- lapply(seq_len(nrow(pts)), function(p) {
      prob <- dbinom(map$s1, n.s[1], pts[p, 1]) * dbinom(map$s2, n.s[2], pts[p, 2])
      data.frame(theta.S1 = pts[p, 1], theta.S2 = pts[p, 2], s1 = map$s1, s2 = map$s2,
                 prob.s = prob, TIE.s = cr$tie.s[, p], contribution = prob * cr$tie.s[, p])
    })
    by.s <- do.call(rbind, parts)
    rownames(by.s) <- NULL
    attr(out, 'by.s') <- by.s
    attr(out, 'TIE') <- data.frame(
      theta.S1 = pts[, 1], theta.S2 = pts[, 2],
      TIE = vapply(parts, function(x) sum(x$contribution), numeric(1))
    )
  }
  if (full) {
    f <- data.frame(s1 = cr$s1, s2 = cr$s2, w1 = cr$w1, w2 = cr$w2, CRP = cr$crp,
                    CRP.total = cr$total)
    f <- f[order(f$s2, f$s1, f$w2, f$w1), , drop = FALSE]
    rownames(f) <- NULL
    attr(out, 'full') <- f
  }
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
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'ss.strata') <- ss.strata
  attr(out, 'search') <- search
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  attr(out, 'reestimation') <- map
  class(out) <- c('bbssr_scrp', 'data.frame')
  out
}
