#' Evaluate One Design of a Grid
#'
#' Internal helper of \code{\link{BinaryGridBSSR}} that evaluates a single design given
#' as a named list of arguments.
#'
#' @param a Named list of arguments of the design
#' @param p Vector of true pooled response probabilities
#' @param Delta.T True treatment effect, or \code{NULL} for the value in \code{a} or the
#'   assumed effect
#' @param type1 Logical. Whether the largest type I error rates are computed
#'
#' @return A list with the initial sample sizes \code{N1} and \code{N2}, the true effect
#'   \code{Delta.T}, the result \code{power} of \code{BinaryPowerBSSR}, its
#'   \code{summary}, and with \code{type1 = TRUE} a one-row data frame \code{type1}
#'
#' @keywords internal
#' @noRd
grid_design <- function(a, p, Delta.T, type1) {
  # Elements are taken with [[ ]], which matches names exactly: $ would match r
  # partially to restricted, rounding or ref.pvalue when r is absent
  for (nm in c('Delta.A', 'r', 'alpha', 'tar.power', 'Test')) {
    if (is.null(a[[nm]])) stop(nm, ' must be given')
  }
  if (xor(is.null(a[['n1.interim']]), is.null(a[['n2.interim']]))) {
    stop('give both n1.interim and n2.interim, or neither')
  }
  if (!is.null(a[['n1.interim']])) {
    a[['n.interim']] <- c(a[['n1.interim']], a[['n2.interim']])
  }
  a[['n1.interim']] <- NULL
  a[['n2.interim']] <- NULL
  if (is.null(Delta.T)) {
    Delta.T <- if (!is.null(a[['Delta.T']])) a[['Delta.T']] else a[['Delta.A']]
  }
  if (length(Delta.T) != 1 || !is.numeric(Delta.T) || is.na(Delta.T)) {
    stop('Delta.T must be a single number')
  }
  a[['Delta.T']] <- NULL
  # Initial sample sizes, given or planned
  if (xor(is.null(a[['N1']]), is.null(a[['N2']]))) stop('give both N1 and N2, or neither')
  if (!is.null(a[['N1']]) && !is.null(a[['p.plan']])) {
    stop('give either N1 and N2 or p.plan, not both')
  }
  if (is.null(a[['N1']])) {
    if (is.null(a[['p.plan']])) stop('give either N1 and N2 or p.plan')
    effect <- if (is.null(a[['effect']])) 'RD' else a[['effect']]
    if (!effect %in% c('RD', 'RR', 'OR')) stop("effect must be 'RD', 'RR' or 'OR'")
    sp <- split_pooled(a[['p.plan']], a[['Delta.A']], a[['r']], effect)
    ss.args <- list(p1 = sp$p1, p2 = sp$p2, r = a[['r']], alpha = a[['alpha']],
                    tar.power = a[['tar.power']], Test = a[['Test']])
    for (nm in c('alternative', 'tsmethod', 'n.grid', 'bb.gamma', 'rounding', 'margin',
                 'ref.pvalue')) {
      if (!is.null(a[[nm]])) ss.args[[nm]] <- a[[nm]]
    }
    if (!is.null(a[['ss.method']])) ss.args[['method']] <- a[['ss.method']]
    ss <- do.call('BinarySampleSize', ss.args)
    a[['N1']] <- ss$N1
    a[['N2']] <- ss$N2
  }
  a[['p.plan']] <- NULL
  power <- do.call('BinaryPowerBSSR', c(list(p = p, Delta.T = Delta.T), a))
  out <- list(N1 = a[['N1']], N2 = a[['N2']], Delta.T = Delta.T, power = power,
              summary = summary(power))
  if (type1) {
    keep <- intersect(names(a), names(formals(BinaryTypeIErrorBSSR)))
    tie <- do.call('BinaryTypeIErrorBSSR', a[keep])
    m <- attr(tie, 'max')
    out[['type1']] <- data.frame(TIE.BSSR = m$TIE[1], theta.BSSR = m$theta[1],
                                 TIE.TRAD = m$TIE[2], theta.TRAD = m$theta[2])
  }
  out
}
