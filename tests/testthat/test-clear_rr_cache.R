test_that("clear_rr_cache removes every stored region", {
  get_rr(5, 4, 0.05, 'Chisq', 'greater', 'minlike', 100, 0)
  expect_gt(.bbssr_cache$cells, 0)
  expect_gt(length(ls(.bbssr_cache$rr)), 0L)
  expect_null(clear_rr_cache())
  expect_equal(.bbssr_cache$cells, 0)
  expect_equal(length(ls(.bbssr_cache$rr)), 0L)
})
