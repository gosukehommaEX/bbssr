test_that("summary.bbssr_grid gives the range of the power and the largest E.N", {
  design <- data.frame(Test = c('Chisq', 'Fisher'))
  res <- BinaryGridBSSR(design, p = c(0.25, 0.4, 0.55), Delta.A = 0.3, N1 = 12, N2 = 12,
                        omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                        ss.method = 'standard')
  s <- summary(res)
  expect_equal(nrow(s), 2)
  expect_named(s, c('design', 'Test', 'N1', 'N2', 'n1.interim', 'n2.interim', 'Delta.T',
                    'power.BSSR.min', 'power.BSSR.max', 'power.TRAD.min',
                    'power.TRAD.max', 'E.N.max'))
  for (i in 1:2) {
    k <- res$design == i
    expect_equal(sum(k), 3)
    expect_identical(s$power.BSSR.min[i], min(res$power.BSSR[k]))
    expect_identical(s$power.BSSR.max[i], max(res$power.BSSR[k]))
    expect_identical(s$power.TRAD.min[i], min(res$power.TRAD[k]))
    expect_identical(s$power.TRAD.max[i], max(res$power.TRAD[k]))
    expect_identical(s$E.N.max[i], max(res$E.N[k]))
  }
  expect_error(summary.bbssr_grid(data.frame(design = 1)), 'no table of designs')
})
