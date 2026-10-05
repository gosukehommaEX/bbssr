# Inputs of max_tail_prob() and max_tail_prob_refined() for the whole outcome grid on a
# given grid of the nuisance parameter, prepared as in unconditional_pvalue()
tail_inputs <- function(stat, N1, N2, theta, decreasing) {
  ord <- order(c(stat), decreasing = decreasing)
  grp <- tie_groups(c(stat)[ord])
  list(
    dbinom1 = outer(0:N1, theta, function(x, t) dbinom(x, N1, t)),
    dbinom2 = outer(0:N2, theta, function(x, t) dbinom(x, N2, t)),
    x1 = as.integer(c(row(stat))[ord] - 1L),
    x2 = as.integer(c(col(stat))[ord] - 1L),
    idx_last = as.integer(grp[, 'last'] - 1L),
    g_lo = rep(0L, length(ord)),
    g_hi = rep(length(theta) - 1L, length(ord)),
    theta = theta, ord = ord, dim = dim(stat)
  )
}

# Tail probabilities in the layout of the outcome grid
tail_max <- function(inp, refined) {
  p <- if (refined) {
    max_tail_prob_refined(inp$dbinom1, inp$dbinom2, inp$x1, inp$x2, inp$idx_last,
                          inp$g_lo, inp$g_hi, inp$theta)
  } else {
    max_tail_prob(inp$dbinom1, inp$dbinom2, inp$x1, inp$x2, inp$idx_last, inp$g_lo,
                  inp$g_hi)
  }
  out <- numeric(length(p))
  out[inp$ord] <- p
  matrix(out, nrow = inp$dim[1], ncol = inp$dim[2])
}

# The expected values in this file are maxima over the unit interval certified by a
# branch and bound on the Bernstein coefficients of the tail probability, computed by
# tools/reference/reference_values.py with numpy and scipy

test_that("refinement never lowers the maximum over the grid", {
  theta <- seq(0, 1, length.out = 100)
  for (tst in c('Z-pool', 'Boschloo')) {
    stat <- if (tst == 'Z-pool') {
      zstat(32, 32)
    } else {
      fisher_pvalue(32, 32, 'greater', 'minlike', FALSE)
    }
    inp <- tail_inputs(stat, 32, 32, theta, decreasing = (tst == 'Z-pool'))
    grid <- tail_max(inp, refined = FALSE)
    ref <- tail_max(inp, refined = TRUE)
    expect_true(all(ref >= grid), info = tst)
    expect_gt(max(ref - grid), 1e-4)
  }
})

test_that("the refined tail probabilities attain the certified maximum", {
  theta <- seq(0, 1, length.out = 100)
  want.sum <- c('Z-pool' = 599.659595259547, Boschloo = 557.854443840286)
  want.flip <- c('Z-pool' = 0.0250058337172828, Boschloo = 0.0250047740156628)
  for (tst in c('Z-pool', 'Boschloo')) {
    stat <- if (tst == 'Z-pool') {
      zstat(32, 32)
    } else {
      fisher_pvalue(32, 32, 'greater', 'minlike', FALSE)
    }
    inp <- tail_inputs(stat, 32, 32, theta, decreasing = (tst == 'Z-pool'))
    grid <- tail_max(inp, refined = FALSE)
    # pmin() takes the attributes of its first argument, so the matrix comes first
    ref <- pmin(tail_max(inp, refined = TRUE), 1)
    expect_equal(sum(ref), want.sum[[tst]], tolerance = 1e-10, info = tst)
    # The outcomes 18 versus 10 and 22 versus 14 are rejected at the level 0.025 on the
    # grid, but their exact p-value exceeds the level
    flip <- cbind(c(18, 22), c(10, 14)) + 1
    expect_true(all(grid[flip] %<<% 0.025), info = tst)
    expect_equal(ref[flip], rep(want.flip[[tst]], 2), tolerance = 1e-10, info = tst)
  }
})

test_that("the refined Berger-Boos p-values attain the certified maximum", {
  stat <- fisher_pvalue(20, 15, 'greater', 'minlike', FALSE)
  p <- unconditional_pvalue(stat, 20, 15, 100, 0.001, decreasing = FALSE,
                            ref.pvalue = TRUE)
  expect_equal(sum(p), 173.890983925883, tolerance = 1e-10)
  grid <- unconditional_pvalue(stat, 20, 15, 100, 0.001, decreasing = FALSE,
                               ref.pvalue = FALSE)
  expect_true(all(p >= grid))
})

test_that("the arcsine grid resolves a maximum missed by the uniform grid", {
  stat <- zstat(150, 60)
  want <- 0.0304116803087155
  # On the uniform grid the largest local maximum of the outcome 103 versus 32 lies
  # between two grid points that are both lower than a neighbouring grid point
  inp <- tail_inputs(stat, 150, 60, seq(0, 1, length.out = 100), decreasing = TRUE)
  expect_lt(tail_max(inp, refined = TRUE)[104, 33], want - 1e-6)
  p <- unconditional_pvalue(stat, 150, 60, 100, 0, decreasing = TRUE, ref.pvalue = TRUE)
  expect_equal(p[104, 33], want, tolerance = 1e-10)
})

test_that("max_tail_prob_refined rejects inconsistent inputs", {
  inp <- tail_inputs(zstat(4, 3), 4, 3, seq(0, 1, length.out = 11), decreasing = TRUE)
  expect_error(max_tail_prob_refined(inp$dbinom1, inp$dbinom2, inp$x1, inp$x2,
                                     inp$idx_last, inp$g_lo, inp$g_hi, inp$theta[-1]),
               'grid does not match')
  expect_error(max_tail_prob_refined(inp$dbinom1, inp$dbinom2, inp$x1, inp$x2[-1],
                                     inp$idx_last, inp$g_lo, inp$g_hi, inp$theta),
               'differ in length')
  bad <- inp$idx_last
  bad[length(bad)] <- 0L
  expect_error(max_tail_prob_refined(inp$dbinom1, inp$dbinom2, inp$x1, inp$x2, bad,
                                     inp$g_lo, inp$g_hi, inp$theta), 'invalid cell index')
  expect_error(max_tail_prob_refined(inp$dbinom1, inp$dbinom2, inp$x1, inp$x2,
                                     inp$idx_last, inp$g_hi, inp$g_lo, inp$theta),
               'invalid grid range')
})
