test_that("print.bbssr_sbssr prints the re-estimated sizes of the strata", {
  re <- StratBSSR(n1 = c(10, 12), n2 = c(10, 12), S = c(9, 16), Delta.A = 0.2, r = 1,
                  alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(20, 24),
                  N2 = c(20, 24), ss.method = 'standard')
  out <- capture.output(r <- print(re))
  expect_identical(r, re)
  expect_true(any(grepl('stratified trial', out)))
  expect_true(any(grepl(sprintf('Final total      : %d', attr(re, 'N.final')), out,
                        fixed = TRUE)))
})

test_that("print.bbssr_sbssr prints a subset without its attributes as a data frame", {
  re <- StratBSSR(n1 = c(10, 12), n2 = c(10, 12), S = c(9, 16), Delta.A = 0.2, r = 1,
                  alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(20, 24),
                  N2 = c(20, 24), ss.method = 'standard')
  # Selecting columns removes the attributes; selecting rows keeps them
  sub <- re[, names(re)]
  expect_null(attr(sub, 'ss.method'))
  out <- capture.output(print(sub))
  expect_false(any(grepl('stratified trial', out)))
  expect_true(any(grepl('N1.final', out, fixed = TRUE)))
})

test_that("print.bbssr_sbssr marks the re-estimation from the pooled estimate", {
  re <- StratBSSR(n1 = c(10, 12), n2 = c(10, 12), S = c(9, 16), Delta.A = 0.2, r = 1,
                  alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(20, 24),
                  N2 = c(20, 24), ss.method = 'standard', ss.strata = 'pooled')
  out <- capture.output(print(re))
  expect_true(any(grepl('(unrestricted, strata pooled)', out, fixed = TRUE)))
  # The pooled estimate is the same in both strata
  expect_equal(re$hat.p, rep(25 / 44, 2))
})
