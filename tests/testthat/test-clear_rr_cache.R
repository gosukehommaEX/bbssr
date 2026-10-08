test_that("clear_rr_cache removes every stored matrix", {
  get_pvalue(5L, 4L, 'Chisq', 'greater', 'minlike', 100L, 0, FALSE, 0, 'RD')
  expect_gt(.bbssr_cache$cells, 0)
  expect_gt(length(ls(.bbssr_cache$pv)), 0L)
  expect_null(clear_rr_cache())
  expect_equal(.bbssr_cache$cells, 0)
  expect_equal(length(ls(.bbssr_cache$pv)), 0L)
})
