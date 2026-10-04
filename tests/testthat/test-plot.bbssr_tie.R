test_that("plot.bbssr_tie returns a ggplot object", {
  res <- BinaryTypeIErrorBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                              alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                              ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.1))
  p <- plot(res)
  expect_s3_class(p, 'ggplot')
  expect_no_error(ggplot2::ggplot_build(p))
  p2 <- plot(res, main = NA, ref.line = NA, colours = c('black', 'grey50'))
  expect_no_error(ggplot2::ggplot_build(p2))
})
