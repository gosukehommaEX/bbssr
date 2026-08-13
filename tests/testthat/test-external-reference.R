# Comparisons against independent implementations. The unconditional tests search over a
# grid of the nuisance parameter, and different packages use different grids, so the
# comparisons use a fine grid on both sides and a loose tolerance.

test_that("the Boschloo p-value agrees with the Exact package", {
  skip_on_cran()
  skip_if_not_installed('Exact')
  N1 <- 10
  N2 <- 10
  p <- attr(BinaryRR(N1, N2, 0.025, 'Boschloo', n.grid = 1000), 'p.value')
  cells <- list(c(8, 2), c(7, 3), c(9, 1), c(6, 4))
  for (cell in cells) {
    x1 <- cell[1]
    x2 <- cell[2]
    tab <- matrix(c(x1, x2, N1 - x1, N2 - x2), nrow = 2)
    ref <- Exact::exact.test(tab, alternative = 'greater', method = 'boschloo',
                             npNumbers = 1000, to.plot = FALSE)$p.value
    expect_equal(p[x1 + 1, x2 + 1], ref, tolerance = 2e-3,
                 info = sprintf('x1 = %d, x2 = %d', x1, x2))
  }
})

test_that("the Z-pooled p-value agrees with the Exact package", {
  skip_on_cran()
  skip_if_not_installed('Exact')
  N1 <- 10
  N2 <- 10
  p <- attr(BinaryRR(N1, N2, 0.025, 'Z-pool', n.grid = 1000), 'p.value')
  cells <- list(c(8, 2), c(7, 3), c(9, 1))
  for (cell in cells) {
    x1 <- cell[1]
    x2 <- cell[2]
    tab <- matrix(c(x1, x2, N1 - x1, N2 - x2), nrow = 2)
    ref <- Exact::exact.test(tab, alternative = 'greater', method = 'z-pooled',
                             npNumbers = 1000, to.plot = FALSE)$p.value
    expect_equal(p[x1 + 1, x2 + 1], ref, tolerance = 2e-3,
                 info = sprintf('x1 = %d, x2 = %d', x1, x2))
  }
})

test_that("the two-sided Boschloo p-value agrees with exact2x2", {
  skip_on_cran()
  skip_if_not_installed('exact2x2')
  N1 <- 10
  N2 <- 10
  p <- attr(BinaryRR(N1, N2, 0.05, 'Boschloo', alternative = 'two.sided',
                     tsmethod = 'central', n.grid = 1000), 'p.value')
  for (cell in list(c(8, 2), c(7, 3))) {
    x1 <- cell[1]
    x2 <- cell[2]
    ref <- exact2x2::boschloo(x1, N1, x2, N2, alternative = 'two.sided',
                              tsmethod = 'central')$p.value
    expect_equal(p[x1 + 1, x2 + 1], ref, tolerance = 2e-3,
                 info = sprintf('x1 = %d, x2 = %d', x1, x2))
  }
})

test_that("the sample size agrees with a direct search over the power function", {
  skip_on_cran()
  for (tst in c('Fisher', 'Boschloo')) {
    ss <- BinarySampleSize(0.6, 0.2, 1, 0.025, 0.8, tst, n.grid = 100)
    powers <- vapply(seq_len(ss$N2 + 3), function(n) {
      BinaryPower(0.6, 0.2, n, n, 0.025, tst, n.grid = 100)$Power
    }, numeric(1))
    expect_gte(powers[ss$N2], 0.8)
    expect_lt(powers[ss$N2 - 1], 0.8)
  }
})
