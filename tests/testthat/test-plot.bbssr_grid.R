test_that("plot.bbssr_grid returns ggplot objects", {
  design <- expand.grid(Test = c('Chisq', 'Fisher'), omega = c(0.3, 0.5),
                        stringsAsFactors = FALSE)
  res <- BinaryGridBSSR(design, p = c(0.3, 0.45), Delta.A = 0.3, N1 = 12, N2 = 12, r = 1,
                        alpha = 0.025, tar.power = 0.8, ss.method = 'standard')
  for (w in c('power', 'E.N')) {
    p <- plot(res, what = w, colour.by = 'Test', facet.by = 'omega')
    expect_s3_class(p, 'ggplot')
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_no_error(ggplot2::ggplot_build(plot(res)))
  expect_no_error(ggplot2::ggplot_build(plot(res, main = NA, sub = NA, ref.line = NA)))
  expect_error(plot(res, colour.by = 'foo'), 'colour.by must name a column')
  expect_error(plot(res, facet.by = c('Test', 'omega')), 'facet.by must name a column')
  # Panels follow the sorted values of the column
  p <- plot(res, facet.by = 'omega')
  expect_identical(levels(p$data$facet), c('omega = 0.3', 'omega = 0.5'))
})
