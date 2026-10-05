test_that("the power curve uses the refinement of the sample size calculation", {
  ss <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, 'Z-pool', ref.pvalue = TRUE)
  seen <- ref_pvalue_calls(plot(ss, N2.range = 20:24))
  expect_equal(nrow(seen), 5L)
  expect_true(all(seen$ref.pvalue))
  ss <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, 'Z-pool')
  seen <- ref_pvalue_calls(plot(ss, N2.range = 20:24))
  expect_false(any(seen$ref.pvalue))
})
