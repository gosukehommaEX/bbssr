#' Power and Sample Size of Several Stratified Re-estimation Designs
#'
#' Evaluates a set of stratified designs with blinded sample size re-estimation over a
#' common set of scenarios of true pooled response probabilities of the strata, by calling
#' \code{\link{StratPowerBSSR}} for each design, and returns the results as one data frame
#' with a row for each design and each scenario.
#'
#' @param design Data frame with one row per design. Its columns are arguments of
#'   \code{\link{StratPowerBSSR}} other than \code{p} and \code{n.interim}, for example
#'   \code{Test}, \code{Delta.A}, \code{omega} or \code{N.max}. An argument with one
#'   value per stratum is given by one column per stratum, named with the suffixes
#'   \code{.S1}, \code{.S2}, ..., for example \code{N1.S1} and \code{N1.S2}; a column
#'   without the suffix gives a value common to the strata, where the argument allows it
#'   (\code{r} and \code{Delta.T}). The interim group sizes are given by the columns
#'   \code{n1.interim.S1, ...} and \code{n2.interim.S1, ...}, and the planned group sizes
#'   either by \code{N1.S1, ...} and \code{N2.S1, ...} or by the planning probabilities
#'   \code{p.plan.S1, ...}, see Details. A missing value in a column means that the
#'   argument is not given for that design
#' @param p True pooled response probabilities of the strata at which every design is
#'   evaluated: a vector with one element per stratum, or a matrix with one row per
#'   scenario and one column per stratum. Every design must have this number of strata.
#'   Scenarios in which a response probability of a group falls outside the unit interval
#'   are dropped for that design, and every design needs at least one scenario that is
#'   kept
#' @param Delta.T True treatment effect, a single value common to the strata or one value
#'   per stratum. The default of \code{NULL} uses the assumed effect \code{Delta.A} of
#'   each design. It can also be given in \code{design}, but not in both places
#' @param type1 Logical. If \code{TRUE}, the largest type I error rate of each design and
#'   of the fixed design is added, as computed by \code{\link{StratTypeIErrorBSSR}} with
#'   its default grid of \code{theta} and, with two strata, the certified maximum. With
#'   three or more strata the rates are simulated at every point of the grid, which takes
#'   long. Default is \code{FALSE}
#' @param verbose Logical. If \code{TRUE}, a message is issued as each design is
#'   evaluated. Default is \code{FALSE}
#' @param ... Arguments shared by all designs, with the names of the arguments of
#'   \code{\link{StratPowerBSSR}}. An argument with one value per stratum is given here as
#'   a vector, for example \code{N1 = c(20, 20)}, and the interim group sizes as
#'   \code{n1.interim} and \code{n2.interim}. An argument may not be given both here and
#'   in \code{design}
#'
#' @return An object of class \code{bbssr_sgrid}, which inherits from \code{bbssr_grid}: a
#'   data frame with one row for each design and each scenario in which the response
#'   probabilities lie in the unit interval, containing the index \code{design} of the
#'   design, the columns of \code{design} (other than the group sizes and \code{Delta.T}),
#'   and:
#' \describe{
#'   \item{N1.S1, ..., N2.S1, ...}{Planned group sizes}
#'   \item{n1.interim.S1, ..., n2.interim.S1, ...}{Interim group sizes}
#'   \item{Delta.T.S1, ...}{True treatment effect in the strata}
#'   \item{scenario}{Row of \code{p}}
#'   \item{p.S1, ..., p1.S1, ..., p2.S1, ...}{True pooled response probabilities of the
#'     strata and those of the two groups}
#'   \item{power.BSSR, power.TRAD}{Power of the design with re-estimation and of the
#'     fixed design}
#'   \item{se.BSSR, se.TRAD}{Monte Carlo standard errors, present when a design has them
#'     (three or more strata, or a random allocation) and \code{NA} for the other designs}
#'   \item{E.N, SD.N}{Expected value and standard deviation of the final total sample
#'     size}
#'   \item{N.q25, N.q50, N.q75}{Quartiles of the final total sample size}
#'   \item{P.N.max}{Probability that the final total sample size is the largest the
#'     design can reach, see \code{\link{summary.bbssr_spowerbssr}}}
#' }
#' The attribute \code{designs} has one row per design with its index, the columns of
#' \code{design} other than the group sizes and \code{Delta.T}, the planned and interim
#' group sizes and the true effect and, with \code{type1 = TRUE}, the largest type I error
#' rates \code{TIE.BSSR} and \code{TIE.TRAD}, the probabilities
#' \code{theta.S1.BSSR, ...} and \code{theta.S1.TRAD, ...} at which they occur, the
#' certified upper bounds \code{bound.BSSR} and \code{bound.TRAD} (\code{NA} without the
#' certificate) and, with a random allocation, the standard errors \code{se.BSSR} and
#' \code{se.TRAD} at the largest rates. Warnings and errors raised while a design is
#' evaluated are prefixed with the index of the design. \code{\link{summary.bbssr_grid}}
#' summarizes each design over the scenarios, and \code{\link{plot.bbssr_sgrid}} draws the
#' results.
#'
#' @details
#' When a design has the planning probabilities \code{p.plan.S1, ...} instead of the
#' planned group sizes, the planning probabilities of the two groups of each stratum are
#' obtained from them and the assumed effect \code{Delta.A}, as in the re-estimation, and
#' the planned group sizes are those of \code{\link{StratSampleSize}} with the fractions
#' \code{frac}, which must then be given, the allocation ratio, the test, the level and
#' the target power of the design, the method \code{ss.method} of the re-estimation and
#' the fixed allocation.
#'
#' Designs with a random allocation (\code{alloc}) or with three or more strata use the
#' random numbers of R in the order of the designs, so a call to \code{set.seed()} before
#' \code{StratGridBSSR} makes the results reproducible.
#'
#' @examples
#' design <- data.frame(Test = c('CMH', 'Strat-Fisher'))
#' grid <- StratGridBSSR(design, p = rbind(c(0.3, 0.5), c(0.4, 0.6)), Delta.A = 0.3,
#'                       N1 = c(10, 10), N2 = c(10, 10), omega = 0.5, r = 1,
#'                       alpha = 0.05, tar.power = 0.8, ss.method = 'standard')
#' summary(grid)
#'
#' \donttest{
#' plot(grid, colour.by = 'Test')
#'
#' # Planned group sizes from the planning probabilities of the strata
#' design <- data.frame(Delta.A = c(0.25, 0.3), p.plan.S1 = 0.3, p.plan.S2 = 0.5)
#' grid <- StratGridBSSR(design, p = c(0.3, 0.5), frac = c(0.5, 0.5), omega = 0.5, r = 1,
#'                       alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                       ss.method = 'standard')
#' attr(grid, 'designs')
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratPowerBSSR}}, \code{\link{StratTypeIErrorBSSR}},
#'   \code{\link{StratSampleSize}}, \code{\link{BinaryGridBSSR}}
#' @export
StratGridBSSR <- function(design, p, Delta.T = NULL, type1 = FALSE, verbose = FALSE, ...) {
  if (is.list(design)) design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (!is.data.frame(design) || nrow(design) < 1) {
    stop('design must be a data frame with at least one row')
  }
  if (!is.numeric(p) || length(p) < 1 || anyNA(p)) {
    stop('p must be a numeric vector or matrix without missing values')
  }
  J <- if (is.matrix(p)) ncol(p) else length(p)
  if (J < 2) stop('p must have one value per stratum, for two or more strata')
  common <- list(...)
  if (length(common) > 0 && (is.null(names(common)) || any(names(common) == ''))) {
    stop('the arguments shared by all designs must be named')
  }
  # Arguments with one value per stratum, given in design by the columns <name>.S1, ...
  per <- c('N1', 'N2', 'n1.interim', 'n2.interim', 'frac', 'r', 'Delta.T', 'p.plan')
  pat <- paste0('^(', gsub('.', '\\.', paste(per, collapse = '|'), fixed = TRUE),
                ')\\.S([0-9]+)$')
  strat.cols <- grep(pat, names(design), value = TRUE)
  base <- sub(pat, '\\1', strat.cols)
  index <- as.integer(sub(pat, '\\2', strat.cols))
  for (b in unique(base)) {
    if (!setequal(index[base == b], seq_len(J)) || anyDuplicated(index[base == b])) {
      stop(sprintf('the columns %s.S1, ... of design must number the strata 1 to %d',
                   b, J))
    }
  }
  given <- c(setdiff(names(design), strat.cols), unique(base))
  allowed <- c(setdiff(names(formals(StratPowerBSSR)), c('p', 'n.interim')),
               'n1.interim', 'n2.interim', 'p.plan')
  unknown <- setdiff(c(given, names(common)), allowed)
  if (length(unknown) > 0) {
    stop('unknown argument(s): ', paste(unknown, collapse = ', '))
  }
  twice <- unique(base[base %in% names(design)])
  if (length(twice) > 0) {
    stop('given both with and without the suffix .S1, ...: ', paste(twice, collapse = ', '))
  }
  both <- intersect(given, names(common))
  if (length(both) > 0) {
    stop('given both in design and as a shared argument: ', paste(both, collapse = ', '))
  }
  if ('Delta.T' %in% given && !is.null(Delta.T)) {
    stop('give Delta.T either as the argument Delta.T or in design')
  }
  sizes <- c('N1', 'N2', 'n1.interim', 'n2.interim', 'Delta.T')
  drop <- c(intersect(names(design), sizes), strat.cols[base %in% sizes])
  design.part <- design[, setdiff(names(design), drop), drop = FALSE]
  value <- function(column, i) {
    v <- column[[i]]
    if (is.factor(v)) as.character(v) else v
  }
  sfx <- function(name) paste0(name, '.S', seq_len(J))
  n.design <- nrow(design)
  rows <- vector('list', n.design)
  info <- vector('list', n.design)
  for (i in seq_len(n.design)) {
    if (verbose) message(sprintf('design %d of %d', i, n.design))
    # Arguments of this design: its row of design, without missing values, with the
    # columns of the strata gathered into vectors, and the shared arguments
    plain <- lapply(design[setdiff(names(design), strat.cols)], value, i = i)
    vals <- plain[!vapply(plain, function(v) length(v) == 1 && is.na(v), logical(1))]
    for (b in unique(base)) {
      v <- vapply(sfx(b), function(cn) as.numeric(design[[cn]][i]), numeric(1))
      if (all(is.na(v))) next
      if (anyNA(v)) {
        stop(sprintf('design %d: the columns %s.S1, ... are missing for some strata',
                     i, b), call. = FALSE)
      }
      vals[[b]] <- unname(v)
    }
    a <- c(vals, common)
    # Errors and warnings are prefixed with the index of the design. A warning issued by
    # the handler is not caught by the same handler
    res <- withCallingHandlers(
      tryCatch(strat_grid_design(a, p, Delta.T, type1), error = function(e) {
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
    sz <- as.data.frame(as.list(c(
      stats::setNames(res[['N1']], sfx('N1')), stats::setNames(res[['N2']], sfx('N2')),
      stats::setNames(attr(pw, 'n1.interim'), sfx('n1.interim')),
      stats::setNames(attr(pw, 'n2.interim'), sfx('n2.interim')),
      stats::setNames(res[['Delta.T']], sfx('Delta.T')))))
    pcols <- c(sfx('p'), sfx('p1'), sfx('p2'))
    se <- data.frame(
      se.BSSR = if (is.null(pw$se.BSSR)) rep(NA_real_, n) else pw$se.BSSR,
      se.TRAD = if (is.null(pw$se.TRAD)) rep(NA_real_, n) else pw$se.TRAD)
    rows[[i]] <- cbind(
      data.frame(design = rep(i, n)),
      design.part[rep(i, n), , drop = FALSE],
      sz[rep(1, n), , drop = FALSE],
      data.frame(scenario = attr(pw, 'scenario')),
      as.data.frame(lapply(pw[pcols], as.numeric)),
      data.frame(power.BSSR = pw$power.BSSR, power.TRAD = pw$power.TRAD), se,
      data.frame(E.N = pw$E.N, SD.N = sm$SD.N, N.q25 = sm$N.q25, N.q50 = sm$N.q50,
                 N.q75 = sm$N.q75, P.N.max = sm$P.N.max)
    )
    info[[i]] <- cbind(data.frame(design = i), design.part[i, , drop = FALSE], sz)
    if (type1) info[[i]] <- cbind(info[[i]], res[['type1']])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  # The standard errors are kept only when some design has them
  for (cn in c('se.BSSR', 'se.TRAD')) {
    if (all(is.na(out[[cn]]))) out[[cn]] <- NULL
  }
  # Only designs with a random allocation have the standard errors of the largest rates;
  # the other designs get NA in them
  cols <- unique(unlist(lapply(info, names)))
  info <- lapply(info, function(d) {
    for (cn in setdiff(cols, names(d))) d[[cn]] <- NA_real_
    d[cols]
  })
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
  class(out) <- c('bbssr_sgrid', 'bbssr_grid', 'data.frame')
  out
}
