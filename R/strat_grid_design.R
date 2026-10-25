#' Evaluate One Stratified Design of a Grid
#'
#' Internal helper of \code{\link{StratGridBSSR}} that evaluates a single design given as
#' a named list of arguments, in which the arguments with one value per stratum are
#' vectors.
#'
#' @param a Named list of arguments of the design
#' @param p Vector or matrix of true pooled response probabilities of the strata
#' @param Delta.T True treatment effect, or \code{NULL} for the value in \code{a} or the
#'   assumed effect
#' @param type1 Logical. Whether the largest type I error rates are computed
#'
#' @return A list with the planned group sizes \code{N1} and \code{N2}, the true effect
#'   \code{Delta.T} with one value per stratum, the result \code{power} of
#'   \code{StratPowerBSSR}, its \code{summary}, and with \code{type1 = TRUE} a one-row
#'   data frame \code{type1}
#'
#' @keywords internal
#' @noRd
strat_grid_design <- function(a, p, Delta.T, type1) {
  # Elements are taken with [[ ]], which matches names exactly: $ would match r
  # partially to restricted when r is absent
  for (nm in c('Delta.A', 'r', 'alpha', 'tar.power', 'Test')) {
    if (is.null(a[[nm]])) stop(nm, ' must be given')
  }
  if (xor(is.null(a[['n1.interim']]), is.null(a[['n2.interim']]))) {
    stop('give both n1.interim and n2.interim, or neither')
  }
  if (!is.null(a[['n1.interim']])) {
    a[['n.interim']] <- rbind(a[['n1.interim']], a[['n2.interim']])
  }
  a[['n1.interim']] <- NULL
  a[['n2.interim']] <- NULL
  if (is.null(Delta.T)) {
    Delta.T <- if (!is.null(a[['Delta.T']])) a[['Delta.T']] else a[['Delta.A']]
  }
  a[['Delta.T']] <- NULL
  # Planned group sizes, given or from the planning probabilities
  if (xor(is.null(a[['N1']]), is.null(a[['N2']]))) stop('give both N1 and N2, or neither')
  if (!is.null(a[['N1']]) && !is.null(a[['p.plan']])) {
    stop('give either N1 and N2 or p.plan, not both')
  }
  if (is.null(a[['N1']])) {
    if (is.null(a[['p.plan']])) stop('give either N1 and N2 or p.plan')
    if (is.null(a[['frac']])) stop('frac must be given with p.plan')
    effect <- if (is.null(a[['effect']])) 'RD' else a[['effect']]
    if (!effect %in% c('RD', 'RR', 'OR')) stop("effect must be 'RD', 'RR' or 'OR'")
    pp <- a[['p.plan']]
    J <- length(pp)
    r <- rep_len(a[['r']], J)
    sp <- lapply(seq_len(J), function(j) split_pooled(pp[j], a[['Delta.A']], r[j], effect))
    ss.args <- list(p1 = vapply(sp, function(s) s$p1, numeric(1)),
                    p2 = vapply(sp, function(s) s$p2, numeric(1)), frac = a[['frac']],
                    r = a[['r']], alpha = a[['alpha']], tar.power = a[['tar.power']],
                    Test = a[['Test']])
    for (nm in c('alternative', 'tsmethod', 'search', 'search.limit')) {
      if (!is.null(a[[nm]])) ss.args[[nm]] <- a[[nm]]
    }
    if (!is.null(a[['ss.method']])) ss.args[['method']] <- a[['ss.method']]
    ss <- do.call('StratSampleSize', ss.args)
    a[['N1']] <- ss$N1
    a[['N2']] <- ss$N2
  }
  a[['p.plan']] <- NULL
  J <- length(a[['N1']])
  if (!is.numeric(Delta.T) || !(length(Delta.T) %in% c(1, J)) || anyNA(Delta.T)) {
    stop('Delta.T must be a single value or one value per stratum')
  }
  power <- do.call('StratPowerBSSR', c(list(p = p, Delta.T = Delta.T), a))
  out <- list(N1 = a[['N1']], N2 = a[['N2']], Delta.T = rep_len(Delta.T, J), power = power,
              summary = summary(power))
  if (type1) {
    keep <- intersect(names(a), names(formals(StratTypeIErrorBSSR)))
    m <- attr(do.call('StratTypeIErrorBSSR', a[keep]), 'max')
    th <- paste0('theta.S', seq_len(J))
    row <- data.frame(TIE.BSSR = m$TIE[1], TIE.TRAD = m$TIE[2])
    for (k in seq_len(J)) {
      row[[paste0(th[k], '.BSSR')]] <- m[[th[k]]][1]
      row[[paste0(th[k], '.TRAD')]] <- m[[th[k]]][2]
    }
    row$bound.BSSR <- if (is.null(m$bound)) NA_real_ else m$bound[1]
    row$bound.TRAD <- if (is.null(m$bound)) NA_real_ else m$bound[2]
    if (!is.null(m$se)) {
      row$se.BSSR <- m$se[1]
      row$se.TRAD <- m$se[2]
    }
    out[['type1']] <- row
  }
  out
}
