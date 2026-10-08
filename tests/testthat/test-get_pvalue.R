test_that("get_pvalue returns the p-values of rr_pvalue", {
  clear_rr_cache()
  for (alt in c('greater', 'less', 'two.sided')) {
    expect_identical(get_pvalue(6L, 4L, 'Z-pool', alt, 'minlike', 30L, 0, FALSE, 0, 'RD'),
                     rr_pvalue(6L, 4L, 'Z-pool', alt, 'minlike', 30L, 0, FALSE, 0, 'RD'),
                     info = alt)
  }
})

test_that("a stored matrix is returned on the second request", {
  clear_rr_cache()
  first <- get_pvalue(9L, 6L, 'Boschloo', 'greater', 'minlike', 50L, 0, FALSE, 0, 'RD')
  expect_equal(.bbssr_cache$cells, 10 * 7)
  second <- get_pvalue(9L, 6L, 'Boschloo', 'greater', 'minlike', 50L, 0, FALSE, 0, 'RD')
  expect_identical(second, first)
  expect_equal(.bbssr_cache$cells, 10 * 7)
  # Any change of an argument gives a separate entry
  get_pvalue(9L, 6L, 'Boschloo', 'greater', 'minlike', 60L, 0, FALSE, 0, 'RD')
  get_pvalue(9L, 6L, 'Boschloo', 'less', 'minlike', 50L, 0, FALSE, 0, 'RD')
  get_pvalue(9L, 6L, 'Boschloo', 'greater', 'minlike', 50L, 0, TRUE, 0, 'RD')
  expect_equal(.bbssr_cache$cells, 4 * 10 * 7)
  expect_equal(length(ls(.bbssr_cache$pv)), 4L)
})

test_that("the store can be switched off", {
  clear_rr_cache()
  old <- options(bbssr.cache = FALSE)
  on.exit(options(old))
  pv <- get_pvalue(6L, 6L, 'Fisher', 'greater', 'minlike', 100L, 0, FALSE, 0, 'RD')
  expect_equal(.bbssr_cache$cells, 0)
  expect_equal(length(ls(.bbssr_cache$pv)), 0L)
  expect_identical(pv, attr(BinaryRR(6, 6, 0.025, 'Fisher'), 'p.value'))
})

test_that("the store is emptied before it exceeds its limit", {
  clear_rr_cache()
  old <- .bbssr_cache$max.cells
  .bbssr_cache$max.cells <- 100
  on.exit({
    .bbssr_cache$max.cells <- old
    clear_rr_cache()
  })
  get_pvalue(6L, 6L, 'Chisq', 'greater', 'minlike', 100L, 0, FALSE, 0, 'RD')
  expect_equal(.bbssr_cache$cells, 49)
  # 49 + 64 cells would exceed the limit, so only the new matrix is kept
  get_pvalue(7L, 7L, 'Chisq', 'greater', 'minlike', 100L, 0, FALSE, 0, 'RD')
  expect_equal(.bbssr_cache$cells, 64)
  expect_equal(length(ls(.bbssr_cache$pv)), 1L)
})

test_that("the scale of the margin is part of the key of a stored matrix", {
  clear_rr_cache()
  rd <- get_pvalue(8L, 7L, 'Farrington-Manning', 'greater', 'minlike', 100L, 0, FALSE,
                   0.5, 'RD')
  rr <- get_pvalue(8L, 7L, 'Farrington-Manning', 'greater', 'minlike', 100L, 0, FALSE,
                   0.5, 'RR')
  expect_equal(length(ls(.bbssr_cache$pv)), 2L)
  expect_false(isTRUE(all.equal(rd, rr)))
  expect_identical(rr, rr_pvalue(8L, 7L, 'Farrington-Manning', 'greater', 'minlike', 100L,
                                 0, FALSE, 0.5, 'RR'))
})
