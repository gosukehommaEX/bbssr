#' Description of a Random Allocation
#'
#' Internal helper returning the text with which the print methods describe the
#' allocation of a two-stratum design whose stratum and group sizes are random.
#'
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param block Size of the permuted blocks, with \code{alloc = 'stratum'}
#' @param n Stratum sizes, with \code{alloc = 'group'}
#'
#' @return A character string
#'
#' @keywords internal
#' @noRd
strat_alloc_label <- function(alloc, block = NULL, n = NULL) {
  switch(alloc,
         stratum = sprintf('random stratum sizes, permuted blocks of %d within the strata',
                           as.integer(block)),
         group = sprintf('stratum sizes (%s), random group sizes',
                         paste(n, collapse = ', ')),
         both = 'random stratum and group sizes')
}
