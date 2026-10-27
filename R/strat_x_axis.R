#' Horizontal Axis of the Plots of Stratified Results
#'
#' Internal helper of \code{\link{plot.bbssr_sgrid}} and
#' \code{\link{plot.bbssr_spowerbssr}} choosing the values on the horizontal axis and its
#' label.
#'
#' @param x Data frame with the columns \code{p.S1, p.S2, ...}
#' @param x.by Name of a numeric column of \code{x}, \code{'scenario'}, or \code{NULL}
#'   for the first of \code{p.S1, p.S2, ...} that varies between the rows, and
#'   \code{'scenario'} when none varies
#' @param xlab Label of the axis, or \code{NULL} for a label made from \code{x.by}
#' @param scenario Index of the scenarios, used for \code{x.by = 'scenario'} when
#'   \code{x} has no column of that name
#'
#' A warning is given when two scenarios of a design (the column \code{design}, if any)
#' share a value on the axis, since the lines of the plot join the scenarios in the order
#' of the axis.
#'
#' @return A list with the values \code{x} on the axis, one per row of \code{x}, and the
#'   label \code{xlab}
#'
#' @keywords internal
#' @noRd
strat_x_axis <- function(x, x.by, xlab, scenario = NULL) {
  if (is.null(x.by)) {
    ps <- grep('^p\\.S[0-9]+$', names(x), value = TRUE)
    varies <- ps[vapply(ps, function(nm) length(unique(x[[nm]])) > 1, logical(1))]
    x.by <- if (length(varies) > 0) varies[1] else 'scenario'
  }
  if (!is.character(x.by) || length(x.by) != 1 || is.na(x.by)) {
    stop('x.by must name a column of x')
  }
  values <- if (x.by %in% names(x)) {
    x[[x.by]]
  } else if (x.by == 'scenario' && !is.null(scenario)) {
    scenario
  } else {
    stop('x.by must name a column of x')
  }
  if (!is.numeric(values)) stop('x.by must name a numeric column of x')
  if (is.null(xlab)) {
    k <- regmatches(x.by, regexec('^p\\.S([0-9]+)$', x.by))[[1]]
    xlab <- if (length(k) == 2) {
      paste('True pooled response probability of stratum', k[2])
    } else if (x.by == 'scenario') {
      'Scenario'
    } else {
      x.by
    }
  }
  design <- if ('design' %in% names(x)) x$design else rep(1L, length(values))
  if (anyDuplicated(data.frame(design, values)) > 0) {
    warning('two scenarios of a design share the value of ', x.by, ' on the horizontal ',
            'axis, so the lines join scenarios that differ in another probability; choose ',
            'x.by, or plot scenarios in which a single probability varies', call. = FALSE)
  }
  list(x = values, xlab = xlab)
}
