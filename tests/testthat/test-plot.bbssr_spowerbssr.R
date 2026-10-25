test_that("plot.bbssr_spowerbssr draws both designs against the varying probability", {
  res <- StratPowerBSSR(p = cbind(0.4, c(0.4, 0.5, 0.6)), Delta.A = 0.3, Delta.T = 0.3,
                        N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05,
                        tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
  g <- plot(res)
  expect_s3_class(g, 'ggplot')
  expect_no_error(ggplot2::ggplot_build(g))
  expect_equal(g$data$p, rep(res$p.S2, 2))
  expect_equal(g$data$Power, c(res$power.BSSR, res$power.TRAD))
  expect_identical(g$labels$x, 'True pooled response probability of stratum 2')
  expect_identical(g$labels$y, 'Power')
  g <- plot(res, x.by = 'scenario', xlab = 'Row of p', ref.line = NA, show.points = FALSE)
  expect_equal(g$data$p, rep(attr(res, 'scenario'), 2))
  expect_identical(g$labels$x, 'Row of p')
  expect_no_error(ggplot2::ggplot_build(g))
  expect_error(plot(res, x.by = 'foo'), 'x.by must name a column')
  # One scenario is drawn at its index
  one <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6),
                        N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                        Test = 'CMH', ss.method = 'standard')
  expect_identical(plot(one)$labels$x, 'Scenario')
})
