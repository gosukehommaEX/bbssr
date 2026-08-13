test_that("tie_groups identifies tie groups of a sorted vector", {
  x <- c(1, 1, 2, 3, 3, 3, 4)
  grp <- tie_groups(x)
  expect_equal(colnames(grp), c('first', 'last'))
  expect_equal(grp[, 'first'], c(1L, 1L, 3L, 4L, 4L, 4L, 7L))
  expect_equal(grp[, 'last'],  c(2L, 2L, 3L, 6L, 6L, 6L, 7L))
})

test_that("tie_groups handles degenerate inputs", {
  expect_equal(nrow(tie_groups(numeric(0))), 0L)
  expect_equal(unname(tie_groups(5)[, 'last']), 1L)
  expect_equal(tie_groups(rep(0, 4))[, 'last'], rep(4L, 4))
})

test_that("tie_groups works on a decreasing sequence", {
  grp <- tie_groups(c(9, 9, 5, 1))
  expect_equal(grp[, 'last'], c(2L, 2L, 3L, 4L))
})

test_that("tie_groups separates values that differ only in relative terms", {
  # An absolute tolerance of 1.5e-08 would merge these two distinct p-values
  grp <- tie_groups(c(1e-12, 2e-12))
  expect_equal(grp[, 'last'], c(1L, 2L))
})

test_that("zstat reproduces the pooled two-sample Z statistic", {
  N1 <- 7
  N2 <- 5
  Z <- zstat(N1, N2)
  expect_equal(dim(Z), c(N1 + 1L, N2 + 1L))
  manual <- outer(0:N1, 0:N2, function(i, j) {
    hat.p <- (i + j) / (N1 + N2)
    se <- sqrt(hat.p * (1 - hat.p) * (1 / N1 + 1 / N2))
    z <- (i / N1 - j / N2) / se
    ifelse(is.finite(z), z, 0)
  })
  expect_equal(Z, manual)
  # Cells with no responders or with all responders have a zero denominator
  expect_equal(Z[1, 1], 0)
  expect_equal(Z[N1 + 1, N2 + 1], 0)
})

test_that("fisher_pvalue matches the hypergeometric tail for a one-sided alternative", {
  N1 <- 6
  N2 <- 8
  p <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  expect_equal(p, fisher_greater_ref(N1, N2))
  expect_false(anyNA(p))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("fisher_pvalue matches stats::fisher.test for a two-sided alternative", {
  N1 <- 6
  N2 <- 7
  p <- fisher_pvalue(N1, N2, 'two.sided', 'minlike', midp = FALSE)
  for (i in 0:N1) {
    for (j in 0:N2) {
      tab <- matrix(c(i, j, N1 - i, N2 - j), nrow = 2)
      expect_equal(p[i + 1, j + 1], stats::fisher.test(tab)$p.value,
                   tolerance = 1e-12, info = sprintf('x1 = %d, x2 = %d', i, j))
    }
  }
})

test_that("fisher_pvalue matches stats::fisher.test for a one-sided alternative", {
  N1 <- 5
  N2 <- 6
  p <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  for (i in 0:N1) {
    for (j in 0:N2) {
      tab <- matrix(c(i, j, N1 - i, N2 - j), nrow = 2)
      expect_equal(p[i + 1, j + 1],
                   stats::fisher.test(tab, alternative = 'greater')$p.value,
                   tolerance = 1e-12)
    }
  }
})

test_that("the central two-sided p-value doubles the smaller tail", {
  N1 <- 5
  N2 <- 5
  p <- fisher_pvalue(N1, N2, 'two.sided', 'central', midp = FALSE)
  for (i in 0:N1) {
    for (j in 0:N2) {
      s <- i + j
      upper <- stats::phyper(i - 1, N1, N2, s, lower.tail = FALSE)
      lower <- stats::phyper(i, N1, N2, s)
      expect_equal(p[i + 1, j + 1], min(1, 2 * min(lower, upper)), tolerance = 1e-12)
    }
  }
})

test_that("the mid-p correction removes half of the observed cell probability", {
  N1 <- 5
  N2 <- 4
  p <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = TRUE)
  manual <- outer(0:N1, 0:N2, function(i, j) {
    stats::phyper(i, N1, N2, i + j, lower.tail = FALSE) +
      0.5 * stats::dhyper(i, N1, N2, i + j)
  })
  expect_equal(p, manual)
})

test_that("the mid-p p-value never exceeds the exact p-value", {
  for (alt in c('greater', 'two.sided')) {
    for (ts in c('minlike', 'central')) {
      exact <- fisher_pvalue(6, 6, alt, ts, midp = FALSE)
      mid <- fisher_pvalue(6, 6, alt, ts, midp = TRUE)
      expect_true(all(mid <= exact + 1e-12))
      expect_true(all(mid >= 0))
    }
  }
})

test_that("cp_bounds returns valid Clopper-Pearson bounds", {
  N <- 20
  gamma <- 0.001
  b <- cp_bounds(N, gamma)
  expect_equal(dim(b), c(N + 1L, 2L))
  expect_equal(unname(b[1, 'lower']), 0)
  expect_equal(unname(b[N + 1, 'upper']), 1)
  expect_true(all(b[, 'lower'] <= b[, 'upper']))
  expect_true(all(diff(b[, 'lower']) >= 0))
  expect_true(all(diff(b[, 'upper']) >= 0))
  # The interval covers the sample proportion
  s <- 0:N
  expect_true(all(b[, 'lower'] <= s / N + 1e-12))
  expect_true(all(b[, 'upper'] >= s / N - 1e-12))
})

test_that("unconditional_pvalue reproduces a direct evaluation of the tail probability", {
  N1 <- 6
  N2 <- 5
  n.grid <- 25
  # Boschloo, ordered by the Fisher p-value with smaller values more extreme
  stat <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  got <- unconditional_pvalue(stat, N1, N2, n.grid, 0, decreasing = FALSE)
  want <- unconditional_ref(stat, N1, N2, n.grid, decreasing = FALSE)
  expect_equal(got, want, tolerance = 1e-12)
  # Z-pooled, ordered by the Z statistic with larger values more extreme
  stat <- zstat(N1, N2)
  got <- unconditional_pvalue(stat, N1, N2, n.grid, 0, decreasing = TRUE)
  want <- unconditional_ref(stat, N1, N2, n.grid, decreasing = TRUE)
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("cells with a tied ordering statistic receive the same p-value", {
  N1 <- 7
  N2 <- 7
  stat <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  p <- unconditional_pvalue(stat, N1, N2, 100, 0, decreasing = FALSE)
  # x1 = 5, x2 = 1 and x1 = 6, x2 = 2 share the Fisher p-value 2 / 39
  expect_equal(stat[6, 2], stat[7, 3], tolerance = 1e-12)
  expect_equal(p[6, 2], p[7, 3], tolerance = 1e-12)
})

test_that("integer_breaks places whole numbers only", {
  br <- integer_breaks()
  expect_equal(br(c(0, 10)), 0:10)
  expect_equal(br(c(-0.5, 10.5)), 0:10)
  # Limits arrive in reverse order when the scale is reversed
  expect_equal(br(c(10.5, -0.5)), 0:10)
  # Long ranges are thinned, and every break is still a whole number
  long <- br(c(0, 200))
  expect_lt(length(long), 21)
  expect_equal(long, round(long))
  expect_true(all(long >= 0 & long <= 200))
  expect_equal(length(br(c(0.2, 0.8))), 0L)
})

test_that("integer_breaks can be clamped to the counts that can occur", {
  # A tile scale reaches half a tile beyond the outermost count and is expanded further,
  # so the limits alone would place breaks at -1 and at N + 1
  br <- integer_breaks(c(0, 10))
  expect_equal(br(c(-1.05, 11.05)), 0:10)
  expect_equal(br(c(11.05, -1.05)), 0:10)
})
