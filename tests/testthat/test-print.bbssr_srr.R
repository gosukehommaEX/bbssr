test_that("print.bbssr_srr prints the settings and the rows with a rejection", {
  rr <- StratRR(c(3, 4), c(2, 5), 0.05, 'Strat-Fisher', 'two.sided', 'central')
  out <- capture.output(res <- print(rr, n = 3))
  expect_identical(res, rr)
  expect_true(any(grepl('Rejection region of a stratified test', out)))
  expect_true(any(grepl('Two-sided method : central', out, fixed = TRUE)))
  expect_true(any(grepl('more rows', out)))
  expect_true(any(grepl('N1 = (3, 4), N2 = (2, 5)', out, fixed = TRUE)))
})
