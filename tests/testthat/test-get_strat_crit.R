test_that("get_strat_crit stores the critical values and returns the same values", {
  clear_rr_cache()
  a <- get_strat_crit(c(5, 6), c(4, 7), 0.05, 'CMH', 'greater', 'blaker')
  expect_named(a, c('N1', 'N2', 'lo', 'hi'))
  # The two-sided method does not matter here, so both calls share one stored result
  b <- get_strat_crit(c(5, 6), c(4, 7), 0.05, 'CMH', 'greater', 'minlike')
  expect_identical(a, b)
  old <- options(bbssr.cache = FALSE)
  on.exit(options(old))
  expect_identical(get_strat_crit(c(5, 6), c(4, 7), 0.05, 'CMH', 'greater', 'minlike'), a)
  crit <- strat_crit(c(5L, 6L), c(4L, 7L), 0L, 0L, 0L, 0.05)
  expect_identical(a$lo, crit$lo)
  expect_identical(a$hi, crit$hi)
})

test_that("get_strat_crit refuses designs with too many stratum totals", {
  expect_error(get_strat_crit(rep(30, 5), rep(30, 5), 0.05, 'CMH', 'greater', 'minlike'),
               'too many')
})
