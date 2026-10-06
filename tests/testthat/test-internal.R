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

test_that("the blaker two-sided p-value follows formula (2) of Mehrotra et al. (2003)", {
  for (n in list(c(6, 9), c(7, 7), c(10, 3))) {
    for (midp in c(FALSE, TRUE)) {
      expect_equal(fisher_pvalue(n[1], n[2], 'two.sided', 'blaker', midp = midp),
                   fisher_blaker_ref(n[1], n[2], midp = midp), tolerance = 1e-12,
                   info = sprintf('%d x %d, midp = %s', n[1], n[2], midp))
    }
  }
})

test_that("the blaker p-value never exceeds the central p-value", {
  # With equal groups the conditional distribution is symmetric and the two coincide, so
  # only unequal groups are used to check that the blaker p-value can be smaller
  for (n in list(c(9, 5), c(14, 7), c(4, 11))) {
    b <- fisher_pvalue(n[1], n[2], 'two.sided', 'blaker', midp = FALSE)
    ce <- fisher_pvalue(n[1], n[2], 'two.sided', 'central', midp = FALSE)
    expect_true(all(b <= ce + 1e-12), info = sprintf('%d x %d', n[1], n[2]))
    expect_gt(sum(b < ce - 1e-6), 0)
  }
})

test_that("the blaker and minlike p-values agree for equal groups only", {
  # The conditional distribution is symmetric when the groups are of equal size, so the
  # two orderings coincide
  for (N in c(6, 11)) {
    expect_equal(fisher_pvalue(N, N, 'two.sided', 'blaker', midp = FALSE),
                 fisher_pvalue(N, N, 'two.sided', 'minlike', midp = FALSE),
                 tolerance = 1e-12, info = N)
  }
  b <- fisher_pvalue(14, 7, 'two.sided', 'blaker', midp = FALSE)
  m <- fisher_pvalue(14, 7, 'two.sided', 'minlike', midp = FALSE)
  expect_gt(max(abs(b - m)), 0.05)
})

test_that("the two-sided Fisher p-values of two published examples are reproduced", {
  # Fay and Hunsberger (2021, Section 8): 8 of 14 against 1 of 7 responders gives 0.087
  # under blaker, 0.159 under minlike and 0.157 under central
  pub <- c(blaker = 0.087, minlike = 0.159, central = 0.157)
  for (ts in names(pub)) {
    p <- fisher_pvalue(14, 7, 'two.sided', ts, midp = FALSE)[9, 2]
    expect_equal(round(p, 3), pub[[ts]], info = ts)
  }
  # The same values to full precision, followed by the mid-p values under blaker and
  # minlike, from tools/reference/reference_values.py
  got <- c(fisher_pvalue(14, 7, 'two.sided', 'blaker', midp = FALSE)[9, 2],
           fisher_pvalue(14, 7, 'two.sided', 'minlike', midp = FALSE)[9, 2],
           fisher_pvalue(14, 7, 'two.sided', 'central', midp = FALSE)[9, 2],
           fisher_pvalue(14, 7, 'two.sided', 'blaker', midp = TRUE)[9, 2],
           fisher_pvalue(14, 7, 'two.sided', 'minlike', midp = TRUE)[9, 2])
  expect_equal(got, c(0.0873065015479876, 0.158823529411765, 0.156656346749226,
                      0.0515479876160991, 0.0873065015479876), tolerance = 1e-10)
  # Mehrotra, Chan and Berger (2003, Section 3.1): 8 of 148 against 1 of 132 gives 0.0388
  p <- fisher_pvalue(148, 132, 'two.sided', 'blaker', midp = FALSE)[9, 2]
  expect_equal(round(p, 4), 0.0388)
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
    for (ts in c('minlike', 'central', 'blaker')) {
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
  got <- unconditional_pvalue(stat, N1, N2, n.grid, 0, decreasing = FALSE,
                              ref.pvalue = FALSE)
  want <- unconditional_ref(stat, N1, N2, n.grid, decreasing = FALSE)
  expect_equal(got, want, tolerance = 1e-12)
  # Z-pooled, ordered by the Z statistic with larger values more extreme
  stat <- zstat(N1, N2)
  got <- unconditional_pvalue(stat, N1, N2, n.grid, 0, decreasing = TRUE,
                              ref.pvalue = FALSE)
  want <- unconditional_ref(stat, N1, N2, n.grid, decreasing = TRUE)
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("cells with a tied ordering statistic receive the same p-value", {
  N1 <- 7
  N2 <- 7
  stat <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  p <- unconditional_pvalue(stat, N1, N2, 100, 0, decreasing = FALSE, ref.pvalue = FALSE)
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
