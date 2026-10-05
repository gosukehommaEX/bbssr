all_tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

test_that("get_rr returns the rejection region of BinaryRR as a plain matrix", {
  clear_rr_cache()
  for (tst in all_tests) {
    for (alt in c('greater', 'less', 'two.sided')) {
      rr <- get_rr(7, 5, 0.05, tst, alt, 'minlike', 30, 0, FALSE)
      expect_identical(rr, as_plain(BinaryRR(7, 5, 0.05, tst, alternative = alt,
                                             n.grid = 30)),
                       info = sprintf('%s, %s', tst, alt))
    }
  }
})

test_that("regions at different levels come from the same stored p-values", {
  clear_rr_cache()
  a <- get_rr(9, 6, 0.025, 'Boschloo', 'greater', 'minlike', 50, 0, FALSE)
  b <- get_rr(9, 6, 0.05, 'Boschloo', 'greater', 'minlike', 50, 0, FALSE)
  expect_equal(length(ls(.bbssr_cache$pv)), 1L)
  # A larger level can only add rejected outcomes
  expect_true(all(b[a]))
  expect_gt(sum(b), sum(a))
})

test_that("the warning on bb.gamma is repeated for a stored region", {
  clear_rr_cache()
  expect_warning(get_rr(5, 5, 0.05, 'Fisher', 'greater', 'minlike', 100, 0.001, FALSE),
                 'bb.gamma')
  expect_warning(get_rr(5, 5, 0.05, 'Fisher', 'greater', 'minlike', 100, 0.001, FALSE),
                 'bb.gamma')
})

test_that("invalid arguments are rejected before anything is stored", {
  clear_rr_cache()
  expect_error(get_rr(c(5, 6), 5, 0.05, 'Chisq', 'greater', 'minlike', 100, 0, FALSE),
               'single value')
  expect_error(get_rr(5.5, 5, 0.05, 'Chisq', 'greater', 'minlike', 100, 0, FALSE),
               'positive integers')
  expect_error(get_rr(5, 5, 1.2, 'Chisq', 'greater', 'minlike', 100, 0, FALSE), 'alpha')
  expect_equal(.bbssr_cache$cells, 0)
})
