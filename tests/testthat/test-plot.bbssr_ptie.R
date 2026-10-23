test_that("plot.bbssr_ptie returns a ggplot object", {
  tie <- PairedTypeIErrorBSSR(Delta.A = 0.25, N = 30, n.interim = 10, alpha = 0.05,
                              tar.power = 0.8, Test = 'McNemar', ss.method = 'standard',
                              psi = seq(0.1, 0.9, by = 0.1), maximize = 'grid')
  p <- plot(tie)
  expect_s3_class(p, 'ggplot')
  expect_no_error(ggplot2::ggplot_build(p))
  p2 <- plot(tie, main = NA, ref.line = NA, colours = c('black', 'grey50'))
  expect_no_error(ggplot2::ggplot_build(p2))
  tie <- PairedTypeIErrorBSSR(Delta.A = 0, N = 30, n.interim = 10, alpha = 0.05,
                              tar.power = 0.8, Test = 'Paired-Score',
                              alternative = 'equivalence', margin = 0.3,
                              ss.method = 'standard', psi = c(0.4, 0.6, 0.8),
                              maximize = 'grid')
  p <- plot(tie)
  expect_s3_class(p$facet, 'FacetWrap')
  expect_no_error(ggplot2::ggplot_build(p))
})
