#' Allocation Sequence of a Two-Stratum Design with Random Sizes
#'
#' Internal helper drawing, with the random numbers of R, the strata and the groups of
#' the patients of a two-stratum trial in the order of enrolment, and returning the sizes
#' of the four groups after every number of patients. The stratum of each patient is drawn
#' with the probability \code{frac[1]} of stratum 1, except with \code{alloc = 'group'},
#' where the first \code{T} patients hold \code{floor(T frac[1])} patients of stratum 1.
#' Within stratum \code{j} the groups follow permuted blocks of size \code{block} with
#' \code{g[j]} patients of group 1 (\code{alloc = 'stratum'}), or each patient is in group
#' 1 with the probability \code{b[j]}. The draws of a stratum do not depend on the length
#' of the sequence, so \code{extend} continues a sequence with the same distribution.
#'
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param al List returned by \code{check_strat_alloc}
#' @param L Number of patients
#' @param x A sequence returned by an earlier call, to be extended to \code{L} patients, or
#'   \code{NULL} (default)
#'
#' @return A list with the strata \code{st} and groups \code{grp} of the patients, the
#'   group sequences \code{gs} of the two strata, and the integer matrix \code{counts}
#'   with one row for each number of patients \code{T = 0, ..., L} and the columns
#'   \code{N1} and \code{N2} of stratum 1 and \code{N1} and \code{N2} of stratum 2
#'
#' @keywords internal
#' @noRd
#' @importFrom stats runif
strat_alloc_seq <- function(alloc, al, L, x = NULL) {
  if (is.null(x)) x <- list(st = integer(0), gs = list(integer(0), integer(0)))
  L0 <- length(x$st)
  if (L > L0) {
    new <- (L0 + 1):L
    x$st <- c(x$st, if (alloc == 'group') {
      # Patients of stratum 1 among the first T patients: floor(T frac[1])
      k <- floor(c(L0, new) * al$frac[1] + 1e-9)
      ifelse(diff(k) == 1, 1L, 2L)
    } else {
      ifelse(runif(length(new)) < al$frac[1], 1L, 2L)
    })
  }
  st <- x$st[seq_len(L)]
  for (j in 1:2) {
    need <- sum(st == j) - length(x$gs[[j]])
    if (need <= 0) next
    x$gs[[j]] <- c(x$gs[[j]], if (alloc == 'stratum') {
      # Whole blocks, each a random permutation of g[j] patients of group 1 and the rest of
      # group 2
      unlist(lapply(seq_len(ceiling(need / al$block)), function(i) {
        sample(rep(1:2, c(al$g[j], al$block - al$g[j])))
      }))
    } else {
      ifelse(runif(need) < al$b[j], 1L, 2L)
    })
  }
  grp <- integer(L)
  for (j in 1:2) {
    idx <- which(st == j)
    grp[idx] <- x$gs[[j]][seq_along(idx)]
  }
  code <- 2L * (st - 1L) + grp
  counts <- if (L == 0) {
    matrix(0L, 1, 4)
  } else {
    rbind(0L, matrix(apply(outer(code, 1:4, '=='), 2, cumsum), ncol = 4))
  }
  storage.mode(counts) <- 'integer'
  dimnames(counts) <- NULL
  x$grp <- grp
  x$counts <- counts
  x
}
