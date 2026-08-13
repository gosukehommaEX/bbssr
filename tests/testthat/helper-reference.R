# Reference implementations used to validate the package internals. They follow the
# definitions of the tests directly and are deliberately slow, so the tests that rely on
# them use very small outcome grids.

# One-sided Fisher exact p-value, computed cell by cell from the hypergeometric tail
fisher_greater_ref <- function(N1, N2) {
  outer(0:N1, 0:N2, function(i, j) stats::phyper(i - 1, N1, N2, i + j, lower.tail = FALSE))
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
