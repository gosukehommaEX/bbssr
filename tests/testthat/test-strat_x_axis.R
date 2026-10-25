test_that("strat_x_axis takes the first probability that varies", {
  x <- data.frame(p.S1 = c(0.4, 0.4, 0.4), p.S2 = c(0.5, 0.6, 0.7),
                  p.S3 = c(0.1, 0.2, 0.3), E.N = c(40, 42, 44))
  ax <- strat_x_axis(x, NULL, NULL)
  expect_identical(ax$x, x$p.S2)
  expect_identical(ax$xlab, 'True pooled response probability of stratum 2')
  ax <- strat_x_axis(x, 'p.S3', 'Stratum 3')
  expect_identical(ax$x, x$p.S3)
  expect_identical(ax$xlab, 'Stratum 3')
  expect_identical(strat_x_axis(x, 'E.N', NULL)$xlab, 'E.N')
  expect_identical(strat_x_axis(x, 'p.S1', NA)$xlab, NA)
})

test_that("strat_x_axis falls back on the scenarios", {
  x <- data.frame(p.S1 = c(0.4, 0.4), p.S2 = c(0.5, 0.5), Test = c('CMH', 'Cochran'))
  ax <- strat_x_axis(x, NULL, NULL, c(2L, 5L))
  expect_identical(ax$x, c(2L, 5L))
  expect_identical(ax$xlab, 'Scenario')
  # A column of that name comes first
  y <- cbind(x, scenario = c(1L, 3L))
  expect_identical(strat_x_axis(y, 'scenario', NULL, c(2L, 5L))$x, c(1L, 3L))
  expect_error(strat_x_axis(x, NULL, NULL), 'x.by must name a column')
  expect_error(strat_x_axis(x, 'foo', NULL), 'x.by must name a column')
  expect_error(strat_x_axis(x, c('p.S1', 'p.S2'), NULL), 'x.by must name a column')
  expect_error(strat_x_axis(x, 'Test', NULL), 'numeric column')
})
