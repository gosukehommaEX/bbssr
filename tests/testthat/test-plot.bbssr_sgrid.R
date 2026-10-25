test_that("plot.bbssr_sgrid returns ggplot objects", {
  design <- data.frame(Test = c('CMH', 'Cochran'), omega = c(0.5, 0.5))
  p <- rbind(c(0.4, 0.5), c(0.4, 0.6), c(0.4, 0.7))
  res <- StratGridBSSR(design, p = p, Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), r = 1,
                       alpha = 0.05, tar.power = 0.8, ss.method = 'standard')
  for (w in c('power', 'E.N')) {
    g <- plot(res, what = w, colour.by = 'Test')
    expect_s3_class(g, 'ggplot')
    expect_no_error(ggplot2::ggplot_build(g))
  }
  # The default axis is the first probability that varies between the scenarios
  g <- plot(res)
  expect_identical(g$data$p[seq_len(nrow(res))], res$p.S2)
  expect_identical(g$labels$x, 'True pooled response probability of stratum 2')
  g <- plot(res, x.by = 'scenario', xlab = 'Row of p')
  expect_identical(g$labels$x, 'Row of p')
  expect_equal(g$data$p[seq_len(nrow(res))], res$scenario)
  expect_no_error(ggplot2::ggplot_build(plot(res, facet.by = 'Test', main = NA)))
  expect_error(plot(res, x.by = 'foo'), 'x.by must name a column')
  expect_error(plot(res, x.by = 'Test'), 'numeric column')
  # A single scenario is drawn against its index
  one <- StratGridBSSR(design, p = c(0.4, 0.6), Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6),
                       r = 1, alpha = 0.05, tar.power = 0.8, ss.method = 'standard')
  expect_identical(plot(one)$labels$x, 'Scenario')
})
