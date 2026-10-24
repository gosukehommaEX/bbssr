#' Final Group Sizes of a Two-Stratum Design for Every Interim Outcome
#'
#' Internal helper applying \code{strat_reestimate} to every pair of interim pooled
#' numbers of responders \code{(s1, s2)} of a design with two strata, and collecting the
#' distinct final designs with their critical values for the compiled code.
#'
#' @param n1 Interim sizes of group 1 in the two strata
#' @param n2 Interim sizes of group 2 in the two strata
#' @param alpha Level of the final analysis
#' @param Test Test of the final analysis
#' @param re.args List of the remaining arguments of \code{strat_reestimate}
#' @param crit Logical. If \code{FALSE}, the designs hold only their group sizes, without
#'   critical values, for callers that need only the final group sizes. Default is
#'   \code{TRUE}
#' @inheritParams StratRR
#'
#' @return A list with the data frame \code{map} (one row per interim outcome, \code{s1}
#'   varying fastest, with the blinded estimates \code{hat.p.S1} and \code{hat.p.S2},
#'   the planning probabilities, \code{n.raw}, \code{n.re}, the final group sizes
#'   \code{N1.S1}, \code{N1.S2}, \code{N2.S1}, \code{N2.S2} and their total \code{N}), the
#'   list \code{designs} of the distinct final designs, with their critical values when
#'   \code{crit = TRUE}, and the
#'   integer vector \code{id} of zero based positions in \code{designs}
#'
#' @keywords internal
#' @noRd
strat_map <- function(n1, n2, alpha, Test, alternative, tsmethod, re.args, crit = TRUE) {
  s <- as.matrix(expand.grid(s1 = 0:(n1[1] + n2[1]), s2 = 0:(n1[2] + n2[2])))
  re <- do.call(strat_reestimate,
                c(list(s = s, n1 = n1, n2 = n2, alternative = alternative,
                       tsmethod = tsmethod), re.args))
  map <- data.frame(s1 = s[, 1], s2 = s[, 2],
                    hat.p.S1 = re$hat.p[, 1], hat.p.S2 = re$hat.p[, 2],
                    hat.p1.S1 = re$hat.p1[, 1], hat.p1.S2 = re$hat.p1[, 2],
                    hat.p2.S1 = re$hat.p2[, 1], hat.p2.S2 = re$hat.p2[, 2],
                    n.raw = re$n.raw, n.re = re$n.re, limit = re$limit,
                    N1.S1 = re$N1[, 1], N1.S2 = re$N1[, 2],
                    N2.S1 = re$N2[, 1], N2.S2 = re$N2[, 2])
  map$N <- as.integer(rowSums(re$N1) + rowSums(re$N2))
  key <- paste(map$N1.S1, map$N1.S2, map$N2.S1, map$N2.S2, sep = '_')
  first <- which(!duplicated(key))
  designs <- lapply(first, function(i) {
    if (!crit) return(list(N1 = re$N1[i, ], N2 = re$N2[i, ]))
    get_strat_crit(re$N1[i, ], re$N2[i, ], alpha, Test, alternative, tsmethod)
  })
  list(map = map, designs = designs, id = as.integer(match(key, key[first]) - 1L))
}
