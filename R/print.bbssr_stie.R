#' Print Method for bbssr_stie Objects
#'
#' Prints the largest type I error rates of a stratified design with blinded sample size
#' re-estimation and of the fixed design, returned by \code{\link{StratTypeIErrorBSSR}}.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_stie}
#' @param digits Number of significant digits. Default is 6
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' tie <- StratTypeIErrorBSSR(Delta.A = 0.2, N1 = c(20, 20), N2 = c(20, 20), omega = 0.5,
#'                            r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                            ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.2))
#' print(tie)
#'
#' @export
print.bbssr_stie <- function(x, digits = 6, ...) {
  if (is.null(attr(x, 'max')) || !all(c('TIE.BSSR', 'TIE.TRAD') %in% names(x))) {
    return(NextMethod())
  }
  cat('Type I error rate of blinded sample size re-estimation for a stratified design\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Planned sizes    : (%s) and (%s)\n', paste(attr(x, 'N1'), collapse = ', '),
              paste(attr(x, 'N2'), collapse = ', ')))
  alloc <- attr(x, 'alloc')
  if (!is.null(alloc) && alloc != 'fixed') {
    cat(sprintf('  Allocation       : %s\n', strat_alloc_label(alloc, attr(x, 'block'))))
    if (is.null(attr(x, 'nseq'))) {
      cat(sprintf('  Interim analysis : first %d patients\n', attr(x, 'n.interim')))
    } else {
      cat(sprintf('  Sequences        : %d, interim analysis of the first %d patients\n',
                  attr(x, 'nseq'), attr(x, 'n.interim')))
    }
  }
  cat(sprintf('  Points           : %d\n', nrow(x)))
  cat(sprintf('  Maximum          : %s\n', maximize_label(x)))
  if (!is.null(attr(x, 'nsim'))) {
    cat(sprintf('  Simulation       : %d trials per point\n', attr(x, 'nsim')))
  }
  cat(sprintf('  Alpha            : %s\n\n', format(attr(x, 'alpha'))))
  mx <- attr(x, 'max')
  tab <- data.frame(Design = ifelse(mx$Design == 'BSSR', 'BSSR', 'Fixed sample'),
                    mx[grep('^theta', names(mx))], TIE = signif(mx$TIE, digits),
                    check.names = FALSE)
  if (any(!is.na(mx$bound))) tab$bound <- signif(mx$bound, digits)
  if (!is.null(mx$se)) tab$se <- signif(mx$se, 2)
  print(tab, row.names = FALSE)
  invisible(x)
}
