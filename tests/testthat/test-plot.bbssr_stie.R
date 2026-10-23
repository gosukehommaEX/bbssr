test_that("plot.bbssr_stie draws the two heat maps", {
  tie <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.2),
                             maximize = 'grid')
  p <- plot(tie)
  expect_s3_class(p, 'ggplot')
  expect_identical(nrow(p$data), 2L * nrow(tie))
  expect_identical(levels(p$data$Design), c('BSSR', 'Fixed sample'))
  bad <- tie
  bad$theta.S3 <- bad$theta.S1
  expect_error(plot(bad), 'two strata')
})
