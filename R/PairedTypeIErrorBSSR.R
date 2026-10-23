#' Type I Error Rate of Blinded Sample Size Re-Estimation for Paired Binary Data
#'
#' Calculates the exact type I error rate of a matched-pairs design with blinded sample
#' size re-estimation, and of the fixed design with the planned number of pairs, over the
#' probability of a discordant pair on the boundary of the null hypothesis, together with
#' the largest rate.
#'
#' @param psi Values of the probability \code{p10 + p01} of a discordant pair at which the
#'   rate is evaluated. Values that do not lie on the boundary of the null hypothesis,
#'   below \code{|margin|}, are dropped. Default is \code{seq(0, 1, by = 0.01)}
#' @param maximize How the largest rate is found. \code{'certified'} (default) takes the
#'   largest value on the grid and certifies the maximum over the interval spanned by the
#'   grid with an upper bound, \code{'refined'} refines the largest values on the grid
#'   between neighbouring grid points, and \code{'grid'} takes the largest value on the
#'   grid. See Details
#' @inheritParams PairedPowerBSSR
#'
#' @return An object of class \code{bbssr_ptie}, a data frame with one row for each value
#'   of \code{psi} on the boundary, containing \code{psi}, the probabilities \code{p10}
#'   and \code{p01} on the boundary, and the rates \code{TIE.BSSR} of the design with
#'   re-estimation and \code{TIE.TRAD} of the fixed design. For a test of equivalence the
#'   rates on the two boundaries are stacked, with the column \code{boundary}. The
#'   attribute \code{max} holds the largest rates \code{TIE} of both designs, the values
#'   \code{psi} at which they occur, the upper bounds \code{bound} of the certified
#'   maximization (otherwise \code{NA}) and, for a test of equivalence, the boundary.
#'   The attribute \code{reestimation} holds the final number of pairs for every interim
#'   number of discordant pairs, and the other attributes the design settings.
#'
#' @details
#' The type I error rate is evaluated on the boundary \code{p10 - p01 = d0} of the null
#' hypothesis, with \code{d0 = 0} for a test of superiority, \code{d0 = -margin} for
#' \code{alternative = 'greater'} and \code{d0 = margin} for \code{'less'}, as for two
#' independent groups. For a fixed design, Berger and Sidik (2003, Theorem 1) show that the
#' largest rate over the null hypothesis of superiority lies on this boundary for rejection
#' regions that are Barnard convex. A test of equivalence has two boundaries,
#' \code{-lower} and \code{upper}, and the largest rate is the larger of the maxima on the
#' two.
#'
#' On a boundary the rejection probability is a polynomial in \code{psi}: given the
#' numbers of discordant pairs, the number of pairs in which only treatment 1 responds is
#' binomial with probability \code{p10 / psi}, which is fixed on the boundary up to the
#' linear reparametrization of \code{paired_boundary}, and the re-estimation depends on the
#' interim data only through the number of discordant pairs. The probability is computed
#' exactly. With \code{maximize = 'certified'} the Bernstein coefficients of the
#' polynomial over the interval spanned by the values of \code{psi} on the boundary are
#' bounded by subdivision, as in \code{\link{BinaryTypeIErrorBSSR}}, which gives an upper
#' bound of the largest rate over the interval within \code{1e-12} of the maximum found.
#'
#' The exact conditional McNemar test keeps its level whatever blinded rule is used
#' (Schepers et al., 2026), so its largest rate does not exceed \code{alpha}.
#'
#' @references
#' Berger RL, Sidik K (2003). Exact unconditional tests for a 2 x 2 matched-pairs design.
#' \emph{Statistical Methods in Medical Research}, 12(2), 91-108.
#'
#' Schepers M, Brannath W, Hoffmann E, Stingl J, Schmidtmann I (2026). Blinded sample size
#' review for McNemar's test based on primary and surrogate endpoints. arXiv:2608.17784.
#'
#' @examples
#' tie <- PairedTypeIErrorBSSR(Delta.A = 0.15, N = 72, n.interim = 30, alpha = 0.05,
#'                             tar.power = 0.8, Test = 'McNemar',
#'                             alternative = 'two.sided', ss.method = 'standard')
#' print(tie)
#'
#' \donttest{
#' plot(tie)
#' }
#'
#' @export
PairedTypeIErrorBSSR <- function(Delta.A, N, omega = NULL, alpha, tar.power, Test,
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
                                 margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  maximize <- match.arg(maximize)
  margin <- norm_margin(alternative, margin, 'RD')
  if (length(psi) < 1 || !is.numeric(psi) || anyNA(psi) || any(psi < 0 | psi > 1)) {
    stop('psi must be a vector of values in [0, 1]')
  }
  psi <- sort(unique(psi))
  if (length(N) != 1 || is.na(N) || N < 1 || N != round(N)) {
    stop('N must be a single positive integer')
  }
  n1 <- paired_n1(N, omega, n.interim)
  # The test of the final analysis, with the Berger-Boos parameter only for the
  # unconditional test
  Test <- check_paired_args(N, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                            margin)$Test
  bb.test <- if (Test == 'McNemar-uncond') bb.gamma else 0
  map <- paired_reestimate(n1, 0:n1, Delta.A, N, tar.power, restricted, alternative,
                           n.grid, bb.gamma, ss.method, ss.Test, ss.alpha, N.min, N.max,
                           ref.pvalue, margin, search, search.limit)
  sides <- null_sides(alternative, margin, 'RD')
  for (s in sides) {
    if (!any(paired_boundary(psi, s$alternative, s$margin)$ok)) {
      stop('no value of psi lies on the boundary of the null hypothesis')
    }
  }
  Ns <- sort(unique(map$N))
  rr.id <- as.integer(match(map$N, Ns) - 1L)
  rr.list <- lapply(Ns, function(n) {
    get_paired_rr(n, alpha, Test, alternative, n.grid, bb.test, ref.pvalue, margin)
  })
  rr.fixed <- get_paired_rr(N, alpha, Test, alternative, n.grid, bb.test, ref.pvalue,
                            margin)
  # Weights of the certified maximization, which do not depend on the boundary. The fixed
  # design has a single final number of pairs and weight 1
  w.bssr <- if (maximize == 'certified') paired_tie_weights(map$N, n1, Ns) else NULL
  per.side <- lapply(sides, function(s) {
    bd <- paired_boundary(psi, s$alternative, s$margin)
    ps <- psi[bd$ok]
    f.bssr <- function(x) {
      b <- paired_boundary(x, s$alternative, s$margin)
      paired_bssr_power(rr.list, rr.id, n1, b$p10, b$p01)
    }
    f.trad <- function(x) {
      b <- paired_boundary(x, s$alternative, s$margin)
      paired_power(rr.fixed, b$p10, b$p01)
    }
    tie.bssr <- f.bssr(ps)
    tie.trad <- f.trad(ps)
    locate <- function(f, y, rr, Nk, W) {
      m <- if (maximize == 'refined') {
        refine_max(f, ps, y)
      } else {
        list(x = ps[which.max(y)], y = max(y))
      }
      m$bound <- NA_real_
      if (maximize == 'certified') {
        cm <- paired_tie_certify(rr, Nk, W, bd$range, s$alternative, s$margin)
        if (cm$y > m$y) {
          m$x <- cm$x
          m$y <- cm$y
        }
        m$bound <- max(cm$bound, m$y)
      }
      m
    }
    b <- paired_boundary(ps, s$alternative, s$margin)
    list(df = data.frame(psi = ps, p10 = b$p10, p01 = b$p01, TIE.BSSR = tie.bssr,
                         TIE.TRAD = tie.trad),
         m.bssr = locate(f.bssr, tie.bssr, rr.list, Ns, w.bssr),
         m.trad = locate(f.trad, tie.trad, list(rr.fixed), N, list(rep(1, N + 1))),
         interval = bd$range, side = s$side)
  })
  # Maxima on the two boundaries that differ by at most 1e-12 count as equal and the lower
  # boundary is reported, as in BinaryTypeIErrorBSSR()
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
  mx <- data.frame(Design = c('BSSR', 'TRAD'), psi = c(m.bssr$x, m.trad$x),
                   TIE = c(m.bssr$y, m.trad$y), bound = c(m.bssr$bound, m.trad$bound),
                   stringsAsFactors = FALSE)
  if (length(per.side) == 1) {
    out <- per.side[[1]]$df
    interval <- per.side[[1]]$interval
  } else {
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
  attr(out, 'N') <- as.integer(N)
  attr(out, 'n.interim') <- n1
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'ss.method') <- ss.method
  attr(out, 'search') <- search
  attr(out, 'maximize') <- maximize
  attr(out, 'nuisance') <- 'psi'
  if (maximize == 'certified') attr(out, 'interval') <- interval
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- ref.pvalue
  attr(out, 'reestimation') <- map
  class(out) <- c('bbssr_ptie', 'data.frame')
  out
}
