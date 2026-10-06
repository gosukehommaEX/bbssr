#' Power and Sample Size of Several Blinded Sample Size Re-estimation Designs
#'
#' Evaluates a set of designs with blinded sample size re-estimation (BSSR) for a binary
#' endpoint over a common set of true pooled response probabilities, by calling
#' \code{\link{BinaryPowerBSSR}} for each design, and returns the results as one data
#' frame with a row for each design and each probability.
#'
#' @param design Data frame with one row per design. Its columns are arguments of
#'   \code{\link{BinaryPowerBSSR}} other than \code{p} and \code{n.interim}, for example
#'   \code{Test}, \code{Delta.A}, \code{omega} or \code{r}. The interim sample sizes are
#'   given by the two columns \code{n1.interim} and \code{n2.interim}, and the initial
#'   sample sizes either by the columns \code{N1} and \code{N2} or by a column
#'   \code{p.plan}, see Details. A missing value in a column means that the argument is
#'   not given for that design, so that designs specified by \code{omega} and by the
#'   interim sample sizes can share one data frame
#' @param p Vector of true pooled response probabilities at which every design is
#'   evaluated. Values at which a response probability falls outside the unit interval are
#'   dropped for that design, and every design needs at least one value that is kept
#' @param Delta.T True treatment effect. The default of \code{NULL} uses the assumed
#'   effect \code{Delta.A} of each design. It can also be given as a column of
#'   \code{design}, but not in both places
#' @param type1 Logical. If \code{TRUE}, the largest type I error rate of each design and
#'   of the corresponding fixed-sample design is added, as computed by
#'   \code{\link{BinaryTypeIErrorBSSR}} with its default grid of \code{theta} and the
#'   certified maximum. It is taken over the common response probability, or over
#'   the boundary of the null hypothesis for a design with a non-inferiority margin.
#'   Default is \code{FALSE}
#' @param verbose Logical. If \code{TRUE}, a message is issued as each design is
#'   evaluated. Default is \code{FALSE}
#' @param ... Arguments shared by all designs, with the same names as the columns of
#'   \code{design}. An argument may not be given both here and as a column
#'
#' @return An object of class \code{bbssr_grid}, a data frame with one row for each
#'   design and each element of \code{p} that gives response probabilities inside the
#'   unit interval, containing the index \code{design} of the design, the columns of
#'   \code{design} (other than the sample sizes and \code{Delta.T}), and:
#' \describe{
#'   \item{N1, N2}{Initial sample sizes}
#'   \item{n1.interim, n2.interim}{Interim sample sizes}
#'   \item{Delta.T}{True treatment effect}
#'   \item{p, p1, p2}{True pooled response probability and those of the two groups}
#'   \item{power.BSSR, power.TRAD}{Power of the BSSR design and of the fixed-sample design
#'     with the initial sample sizes}
#'   \item{E.N, SD.N}{Expected value and standard deviation of the final total sample
#'     size}
#'   \item{N.q25, N.q50, N.q75}{Quartiles of the final total sample size}
#'   \item{P.N.max}{Probability that the final total sample size is the largest the
#'     design can reach}
#' }
#' The attribute \code{designs} has one row per design with its index, the columns of
#' \code{design} other than the sample sizes and \code{Delta.T}, the initial and interim
#' sample sizes and the true effect and, with \code{type1 = TRUE}, the largest type I
#' error rates \code{TIE.BSSR} and \code{TIE.TRAD} and the values \code{theta.BSSR} and
#' \code{theta.TRAD} at which they occur. Warnings and errors raised while a design is
#' evaluated are prefixed with the index of the design.
#'
#' @details
#' When \code{design} has a column \code{p.plan} instead of the initial sample sizes, the
#' planning proportions of the two groups are obtained from the pooled proportion
#' \code{p.plan} and the assumed effect \code{Delta.A}, as in the re-estimation, and the
#' initial sample sizes are those of \code{\link{BinarySampleSize}} with the test, the
#' level, the target power and the allocation ratio of the design, the method
#' \code{ss.method} of the re-estimation and the rounding rule of the design. The test
#' \code{Test} and the level \code{alpha} of the final analysis are used even when the
#' re-estimation uses another test \code{ss.Test} or level \code{ss.alpha}.
#'
#' The quartiles and \code{P.N.max} are those of \code{\link{summary.bbssr_powerbssr}}.
#' The p-values of each test are kept for the rest of the session, so designs that share
#' sample sizes are evaluated faster after the first one.
#'
#' @examples
#' design <- expand.grid(Test = c('Chisq', 'Fisher'), omega = c(0.3, 0.5),
#'                       stringsAsFactors = FALSE)
#' grid <- BinaryGridBSSR(design, p = c(0.3, 0.4, 0.5), Delta.A = 0.3, N1 = 30, N2 = 30,
#'                        r = 1, alpha = 0.025, tar.power = 0.8, ss.method = 'standard')
#' summary(grid)
#'
#' \donttest{
#' plot(grid, colour.by = 'Test', facet.by = 'omega')
#'
#' # Initial sample sizes from a planning proportion, and the largest type I error rates
#' design <- data.frame(Delta.A = c(0.2, 0.3), p.plan = 0.4)
#' grid <- BinaryGridBSSR(design, p = c(0.3, 0.5), omega = 0.5, r = 1, alpha = 0.025,
#'                        tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
#'                        type1 = TRUE)
#' attr(grid, 'designs')
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryPowerBSSR}}, \code{\link{BinaryTypeIErrorBSSR}},
#'   \code{\link{BinarySampleSize}}
#' @export
BinaryGridBSSR <- function(design, p, Delta.T = NULL, type1 = FALSE, verbose = FALSE,
                           ...) {
  if (is.list(design)) design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (!is.data.frame(design) || nrow(design) < 1) {
    stop('design must be a data frame with at least one row')
  }
  if (!is.numeric(p) || length(p) < 1 || anyNA(p)) {
    stop('p must be a numeric vector without missing values')
  }
  common <- list(...)
  allowed <- c(setdiff(names(formals(BinaryPowerBSSR)), c('p', 'n.interim')),
               'n1.interim', 'n2.interim', 'p.plan')
  if (length(common) > 0 && (is.null(names(common)) || any(names(common) == ''))) {
    stop('the arguments shared by all designs must be named')
  }
  unknown <- setdiff(c(names(design), names(common)), allowed)
  if (length(unknown) > 0) {
    stop('unknown argument(s): ', paste(unknown, collapse = ', '))
  }
  both <- intersect(names(design), names(common))
  if (length(both) > 0) {
    stop('given both as a column of design and as a shared argument: ',
         paste(both, collapse = ', '))
  }
  if ('Delta.T' %in% names(design) && !is.null(Delta.T)) {
    stop('give Delta.T either as the argument Delta.T or as a column of design')
  }
  sizes <- c('N1', 'N2', 'n1.interim', 'n2.interim', 'Delta.T')
  design.part <- design[, setdiff(names(design), sizes), drop = FALSE]
  n.design <- nrow(design)
  rows <- vector('list', n.design)
  info <- vector('list', n.design)
  for (i in seq_len(n.design)) {
    if (verbose) message(sprintf('design %d of %d', i, n.design))
    # Arguments of this design: its row of design, without missing values, and the
    # shared arguments
    vals <- lapply(design, function(column) {
      v <- column[[i]]
      if (is.factor(v)) as.character(v) else v
    })
    vals <- vals[!vapply(vals, function(v) length(v) == 1 && is.na(v), logical(1))]
    a <- c(vals, common)
    # Errors and warnings are prefixed with the index of the design. A warning issued by
    # the handler is not caught by the same handler
    res <- withCallingHandlers(
      tryCatch(grid_design(a, p, Delta.T, type1), error = function(e) {
        stop(sprintf('design %d: %s', i, conditionMessage(e)), call. = FALSE)
      }),
      warning = function(w) {
        warning(sprintf('design %d: %s', i, conditionMessage(w)), call. = FALSE)
        invokeRestart('muffleWarning')
      }
    )
    pw <- res[['power']]
    sm <- res[['summary']]
    n <- nrow(pw)
    rows[[i]] <- cbind(
      data.frame(design = rep(i, n)),
      design.part[rep(i, n), , drop = FALSE],
      data.frame(N1 = res[['N1']], N2 = res[['N2']],
                 n1.interim = attr(pw, 'n1.interim'), n2.interim = attr(pw, 'n2.interim'),
                 Delta.T = res[['Delta.T']], p = pw$p, p1 = pw$p1, p2 = pw$p2,
                 power.BSSR = pw$power.BSSR, power.TRAD = pw$power.TRAD, E.N = pw$E.N,
                 SD.N = sm$SD.N, N.q25 = sm$N.q25, N.q50 = sm$N.q50, N.q75 = sm$N.q75,
                 P.N.max = sm$P.N.max)
    )
    info[[i]] <- cbind(
      data.frame(design = i),
      design.part[i, , drop = FALSE],
      data.frame(N1 = res[['N1']], N2 = res[['N2']], n1.interim = attr(pw, 'n1.interim'),
                 n2.interim = attr(pw, 'n2.interim'), Delta.T = res[['Delta.T']])
    )
    if (type1) info[[i]] <- cbind(info[[i]], res[['type1']])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  designs <- do.call(rbind, info)
  rownames(designs) <- NULL
  tar.power <- if ('tar.power' %in% names(design)) {
    design[['tar.power']]
  } else {
    common[['tar.power']]
  }
  attr(out, 'designs') <- designs
  attr(out, 'common') <- common
  attr(out, 'tar.power') <- unique(tar.power)
  attr(out, 'type1') <- type1
  class(out) <- c('bbssr_grid', 'data.frame')
  out
}
