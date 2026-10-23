test_that("paired_tail_max maximizes the tail probabilities over the grid", {
  N <- 8
  Z <- paired_score_z(N, 0)
  ok <- which(!is.na(Z))
  a <- row(Z)[ok] - 1L
  b <- col(Z)[ok] - 1L
  ord <- order(Z[ok], decreasing = TRUE)
  grp <- tie_groups(Z[ok][ord])
  theta <- seq(0, 1, length.out = 51)
  dd <- outer(0:N, theta, function(x, t) stats::dbinom(x, N, t))
  D <- a[ord] + b[ord]
  w <- stats::dbinom(a[ord], D, 0.5)
  p <- paired_tail_max(dd, as.integer(D), w, as.integer(grp[, 'last'] - 1L),
                       rep(0L, length(ord)), rep(50L, length(ord)), theta, FALSE)
  # Direct computation: probability of the outcomes at least as extreme
  for (k in c(1, 5, 20, length(ord))) {
    s <- Z[ok][ord][k]
    tail <- Z[ok] >= s - 1e-10 * pmax(abs(s), abs(Z[ok]))
    prob <- vapply(theta, function(t) {
      sum(stats::dbinom(a + b, N, t)[tail] * stats::dbinom(a, a + b, 0.5)[tail])
    }, numeric(1))
    expect_equal(p[k], max(prob), tolerance = 1e-13)
  }
  # Refinement never lowers the value
  q <- paired_tail_max(dd, as.integer(D), w, as.integer(grp[, 'last'] - 1L),
                       rep(0L, length(ord)), rep(50L, length(ord)), theta, TRUE)
  expect_true(all(q >= p - 1e-15))
})

test_that("paired_tail_max rejects tie groups that would not end", {
  dd <- matrix(0.5, 2, 2)
  expect_error(paired_tail_max(dd, c(0L, 1L), c(1, 1), c(0L, 0L), c(0L, 0L), c(1L, 1L),
                               c(0, 1), FALSE), 'idx_last')
  expect_error(paired_tail_max(dd, c(0L, 2L), c(1, 1), c(0L, 1L), c(0L, 0L), c(1L, 1L),
                               c(0, 1), FALSE), 'D must lie')
})
