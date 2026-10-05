test_that("plot.bbssr_crp returns a ggplot object for each quantity", {
  res <- BinaryCondRejectBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
                              alpha = 0.025, tar.power = 0.8, Test = 'Fisher',
                              ss.method = 'standard', theta = c(0.3, 0.5))
  for (w in c('CRP', 'CRP.total', 'by.s')) {
    p <- plot(res, what = w)
    expect_s3_class(p, 'ggplot')
    expect_no_error(ggplot2::ggplot_build(p))
  }
  p2 <- plot(res, main = NA, sub = NA, colours = c('black', 'grey50'))
  expect_no_error(ggplot2::ggplot_build(p2))
  res0 <- BinaryCondRejectBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6),
                               r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Fisher',
                               ss.method = 'standard')
  expect_error(plot(res0, what = 'by.s'), 'theta')
})
