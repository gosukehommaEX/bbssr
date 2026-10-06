# Reference implementations used to validate the package internals. They follow the
# definitions of the tests directly and are deliberately slow, so the tests that rely on
# them use very small outcome grids.

# One-sided Fisher exact p-value, computed cell by cell from the hypergeometric tail
fisher_greater_ref <- function(N1, N2) {
  outer(0:N1, 0:N2, function(i, j) stats::phyper(i - 1, N1, N2, i + j, lower.tail = FALSE))
}

# Two-sided Fisher p-value of the blaker convention, formula (2) of Mehrotra, Chan and
# Berger (2003), computed table by table from the hypergeometric tails. With midp = TRUE
# the tables tied with the observed one, the observed one included, contribute half of
# their probability
fisher_blaker_ref <- function(N1, N2, midp = FALSE) {
  p <- matrix(NA_real_, nrow = N1 + 1, ncol = N2 + 1)
  for (i in 0:N1) {
    for (j in 0:N2) {
      s <- i + j
      k <- max(0, s - N2):min(N1, s)
      f <- stats::dhyper(k, N1, N2, s)
      g <- pmin(stats::phyper(k, N1, N2, s),
                stats::phyper(k - 1, N1, N2, s, lower.tail = FALSE))
      g0 <- g[k == i]
      tied <- abs(g - g0) <= 1e-10 * pmax(g, g0)
      more <- g < g0 & !tied
      p[i + 1, j + 1] <- min(1, sum(f[more]) + (if (midp) 0.5 else 1) * sum(f[tied]))
    }
  }
  p
}

# Null tail probability of an ordering statistic, maximized over the nuisance parameter,
# obtained by forming the tail set of every cell explicitly
unconditional_ref <- function(stat, N1, N2, n.grid, decreasing) {
  theta <- seq(0, 1, length.out = n.grid)
  joint <- lapply(theta, function(t) {
    outer(stats::dbinom(0:N1, N1, t), stats::dbinom(0:N2, N2, t))
  })
  p <- matrix(0, nrow = N1 + 1, ncol = N2 + 1)
  for (i in 0:N1) {
    for (j in 0:N2) {
      s0 <- stat[i + 1, j + 1]
      tol <- 1e-10 * pmax(abs(stat), abs(s0))
      mask <- if (decreasing) stat >= s0 - tol else stat <= s0 + tol
      p[i + 1, j + 1] <- max(vapply(joint, function(m) sum(m * mask), numeric(1)))
    }
  }
  pmin(p, 1)
}

# Largest null rejection probability of a rejection region over the nuisance parameter
max_type1 <- function(RR, n.grid = 401) {
  N1 <- attr(RR, 'N1')
  N2 <- attr(RR, 'N2')
  m <- matrix(as.vector(RR), nrow = N1 + 1L, ncol = N2 + 1L)
  theta <- seq(0, 1, length.out = n.grid)
  max(vapply(theta, function(t) {
    sum(outer(stats::dbinom(0:N1, N1, t), stats::dbinom(0:N2, N2, t)) * m)
  }, numeric(1)))
}

# Plain logical matrix without the class and attributes of a bbssr_rr object
as_plain <- function(RR) {
  matrix(as.vector(RR), nrow = attr(RR, 'N1') + 1L, ncol = attr(RR, 'N2') + 1L)
}

# Berger-Boos p-value obtained by forming the tail set and the confidence interval of
# every cell explicitly. The grid matches the one built inside unconditional_pvalue.
berger_boos_ref <- function(stat, N1, N2, n.grid, gamma, decreasing) {
  bnd <- cp_bounds(N1 + N2, gamma)
  theta <- sort(unique(c(seq(0, 1, length.out = n.grid), bnd[, 'lower'], bnd[, 'upper'])))
  joint <- lapply(theta, function(t) {
    outer(stats::dbinom(0:N1, N1, t), stats::dbinom(0:N2, N2, t))
  })
  p <- matrix(0, nrow = N1 + 1, ncol = N2 + 1)
  for (i in 0:N1) {
    for (j in 0:N2) {
      s0 <- stat[i + 1, j + 1]
      tol <- 1e-10 * pmax(abs(stat), abs(s0))
      mask <- if (decreasing) stat >= s0 - tol else stat <= s0 + tol
      keep <- which(theta >= bnd[i + j + 1, 'lower'] & theta <= bnd[i + j + 1, 'upper'])
      p[i + 1, j + 1] <- gamma +
        max(vapply(joint[keep], function(m) sum(m * mask), numeric(1)))
    }
  }
  pmin(p, 1)
}

# Rejection probability of a re-estimation design obtained by summing the joint
# probability of the second-stage outcomes over the part of the final rejection region
# that can be reached from each interim cell. This is the summation of version 2.0.0
bssr_power_ref <- function(rr.list, rr.id, x11, x12, n21, n22, p1, p2, n11, n12) {
  vapply(seq_along(p1), function(s) {
    sum(vapply(seq_along(x11), function(cell) {
      k <- rr.id[cell] + 1L
      sub <- rr.list[[k]][x11[cell] + 0:n21[k] + 1L, x12[cell] + 0:n22[k] + 1L,
                          drop = FALSE]
      cp <- sum(outer(stats::dbinom(0:n21[k], n21[k], p1[s]),
                      stats::dbinom(0:n22[k], n22[k], p2[s])) * sub)
      stats::dbinom(x11[cell], n11, p1[s]) * stats::dbinom(x12[cell], n12, p2[s]) * cp
    }, numeric(1)))
  }, numeric(1))
}

# Number of runs of rejected cells in each column of a logical matrix
column_run_count <- function(rr) {
  apply(rr, 2, function(v) sum(diff(c(FALSE, v)) == 1))
}

# Test, refinement flag, margin and two-sided convention of every call of get_pvalue()
# made while code is evaluated, which shows whether ref.pvalue, margin and tsmethod reach
# each p-value matrix used by a function
ref_pvalue_calls <- function(code) {
  seen <- data.frame(Test = character(0), ref.pvalue = logical(0), margin = numeric(0),
                     tsmethod = character(0))
  real <- get_pvalue
  testthat::local_mocked_bindings(
    get_pvalue = function(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma,
                          ref.pvalue, margin) {
      seen[nrow(seen) + 1L, ] <<- list(Test, ref.pvalue, margin, tsmethod)
      real(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma, ref.pvalue, margin)
    }
  )
  force(code)
  seen
}

# Conditional rejection probabilities of a re-estimation design under equal response
# probabilities, summed outcome by outcome over the responder counts of group 1 in the
# two stages, together with the rejection probability given the total alone
bssr_cond_reject_ref <- function(rr.list, rr.id, n21, n22, n11, n12) {
  s.out <- integer(0)
  s2.out <- integer(0)
  crp <- numeric(0)
  crp.total <- numeric(0)
  for (s in 0:(n11 + n12)) {
    k <- rr.id[s + 1L] + 1L
    rr <- rr.list[[k]]
    m1 <- n21[k]
    m2 <- n22[k]
    for (s2 in 0:(m1 + m2)) {
      v <- 0
      for (a in 0:n11) {
        for (b in 0:m1) {
          if (s - a < 0 || s - a > n12 || s2 - b < 0 || s2 - b > m2) next
          if (rr[a + b + 1L, s - a + s2 - b + 1L]) {
            v <- v + stats::dhyper(a, n11, n12, s) * stats::dhyper(b, m1, m2, s2)
          }
        }
      }
      t <- s + s2
      x1 <- max(0, t - n12 - m2):min(n11 + m1, t)
      vt <- sum(stats::dhyper(x1, n11 + m1, n12 + m2, t) *
                  rr[cbind(x1 + 1L, t - x1 + 1L)])
      s.out <- c(s.out, s)
      s2.out <- c(s2.out, s2)
      crp <- c(crp, v)
      crp.total <- c(crp.total, vt)
    }
  }
  list(s = s.out, s2 = s2.out, crp = crp, crp.total = crp.total)
}
