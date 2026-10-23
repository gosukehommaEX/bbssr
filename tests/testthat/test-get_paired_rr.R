test_that("get_paired_rr rejects the p-values below the level", {
  rr <- get_paired_rr(12, 0.05, 'McNemar', 'greater', 100, 0, FALSE, 0)
  pv <- paired_pvalue(12L, 'McNemar', 'greater', 100L, 0, FALSE, 0)
  expect_type(rr, 'logical')
  expect_equal(dim(rr), c(13L, 13L))
  expect_false(anyNA(rr))
  expect_identical(rr[!is.na(pv)], pv[!is.na(pv)] %<<% 0.05)
  # Impossible outcomes are never rejected
  expect_false(any(rr[is.na(pv)]))
  expect_equal(sum(rr), 21)
})

test_that("get_paired_rr validates its arguments", {
  expect_error(get_paired_rr(12, 0.05, 'McNemar', 'greater', 100, 0, FALSE, 0.1),
               'non-zero margin')
  expect_error(get_paired_rr(0, 0.05, 'McNemar', 'greater', 100, 0, FALSE, 0),
               'positive integer')
})
