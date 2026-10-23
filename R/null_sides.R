#' One-Sided Components of the Null Hypothesis
#'
#' Internal helper splitting the null hypothesis into its one-sided components. A test of
#' equivalence has two. The lower one, \code{p1 - p2 <= -lower} (or
#' \code{p1 / p2 <= lower} on the scale of the risk ratio), is the null hypothesis of the
#' one-sided test against \code{'greater'} with the margin \code{lower}, and the upper
#' one, \code{p1 - p2 >= upper} (or \code{p1 / p2 >= upper}), that of the one-sided test
#' against \code{'less'} with the margin \code{upper}. Every other alternative has a
#' single component, the null hypothesis itself. The type I error rate is evaluated on
#' the boundary of each component.
#'
#' @param alternative \code{'greater'}, \code{'less'}, \code{'two.sided'} or
#'   \code{'equivalence'}
#' @param margin Margin, see \code{norm_margin}
#' @param margin.scale \code{'RD'} or \code{'RR'}
#'
#' @return A list of components, each a list with the one-sided \code{alternative}, the
#'   single \code{margin} and the label \code{side}, \code{'lower'}, \code{'upper'} or
#'   \code{NA}
#'
#' @keywords internal
#' @noRd
null_sides <- function(alternative, margin, margin.scale) {
  if (!identical(alternative, 'equivalence')) {
    return(list(list(alternative = alternative, margin = margin, side = NA_character_)))
  }
  m <- norm_margin(alternative, margin, margin.scale)
  list(list(alternative = 'greater', margin = m[1], side = 'lower'),
       list(alternative = 'less', margin = m[2], side = 'upper'))
}
