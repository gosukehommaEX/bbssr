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

test_that("plot.bbssr_tie draws a panel for each boundary of a test of equivalence", {
  tie <- BinaryTypeIErrorBSSR(Delta.A = 0, N1 = 30, N2 = 30, n.interim = c(15, 15),
                              r = 1, alpha = 0.05, tar.power = 0.8,
                              Test = 'Farrington-Manning', alternative = 'equivalence',
                              ss.method = 'standard', margin = 0.3,
                              theta = c(0.3, 0.5, 0.7), maximize = 'grid')
  p <- plot(tie)
  expect_s3_class(p, 'ggplot')
  expect_s3_class(p$facet, 'FacetWrap')
  expect_no_error(ggplot2::ggplot_build(p))
})
