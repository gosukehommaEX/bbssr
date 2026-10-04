all_tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

test_that("get_rr returns the rejection region of BinaryRR as a plain matrix", {
  clear_rr_cache()
  for (tst in all_tests) {
    for (alt in c('greater', 'two.sided')) {
      rr <- get_rr(7, 5, 0.05, tst, alt, 'minlike', 30, 0)
      expect_identical(rr, as_plain(BinaryRR(7, 5, 0.05, tst, alternative = alt,
                                             n.grid = 30)),
                       info = sprintf('%s, %s', tst, alt))
    }
  }
})

test_that("a stored region is returned on the second request", {
  clear_rr_cache()
  first <- get_rr(9, 6, 0.025, 'Boschloo', 'greater', 'minlike', 50, 0)
  expect_equal(.bbssr_cache$cells, 10 * 7)
  second <- get_rr(9, 6, 0.025, 'Boschloo', 'greater', 'minlike', 50, 0)
  expect_identical(second, first)
  expect_equal(.bbssr_cache$cells, 10 * 7)
  # Any change of an argument gives a separate entry
  get_rr(9, 6, 0.025, 'Boschloo', 'greater', 'minlike', 60, 0)
  get_rr(9, 6, 0.05, 'Boschloo', 'greater', 'minlike', 50, 0)
  expect_equal(.bbssr_cache$cells, 3 * 10 * 7)
  expect_equal(length(ls(.bbssr_cache$rr)), 3L)
})

test_that("the store can be switched off", {
  clear_rr_cache()
  old <- options(bbssr.cache = FALSE)
  on.exit(options(old))
  rr <- get_rr(6, 6, 0.025, 'Fisher', 'greater', 'minlike', 100, 0)
  expect_equal(.bbssr_cache$cells, 0)
  expect_equal(length(ls(.bbssr_cache$rr)), 0L)
  expect_identical(rr, as_plain(BinaryRR(6, 6, 0.025, 'Fisher')))
})

test_that("the store is emptied before it exceeds its limit", {
  clear_rr_cache()
  old <- .bbssr_cache$max.cells
  .bbssr_cache$max.cells <- 100
  on.exit({
    .bbssr_cache$max.cells <- old
    clear_rr_cache()
  })
  get_rr(6, 6, 0.025, 'Chisq', 'greater', 'minlike', 100, 0)
  expect_equal(.bbssr_cache$cells, 49)
  # 49 + 64 cells would exceed the limit, so only the new region is kept
  get_rr(7, 7, 0.025, 'Chisq', 'greater', 'minlike', 100, 0)
  expect_equal(.bbssr_cache$cells, 64)
  expect_equal(length(ls(.bbssr_cache$rr)), 1L)
})

test_that("the warning on bb.gamma is repeated for a stored region", {
  clear_rr_cache()
  expect_warning(get_rr(5, 5, 0.05, 'Fisher', 'greater', 'minlike', 100, 0.001),
                 'bb.gamma')
  expect_warning(get_rr(5, 5, 0.05, 'Fisher', 'greater', 'minlike', 100, 0.001),
                 'bb.gamma')
})

test_that("invalid arguments are still rejected by BinaryRR", {
  clear_rr_cache()
  expect_error(get_rr(c(5, 6), 5, 0.05, 'Chisq', 'greater', 'minlike', 100, 0),
               'single value')
  expect_error(get_rr(5.5, 5, 0.05, 'Chisq', 'greater', 'minlike', 100, 0),
               'positive integers')
  expect_equal(.bbssr_cache$cells, 0)
})
