test_that("print.bbssr_srr prints the settings and the rows with a rejection", {
  rr <- StratRR(c(3, 4), c(2, 5), 0.05, 'Strat-Fisher', 'two.sided', 'central')
  out <- capture.output(res <- print(rr, n = 3))
  expect_identical(res, rr)
  expect_true(any(grepl('Rejection region of a stratified test', out)))
  expect_true(any(grepl('Two-sided method : central', out, fixed = TRUE)))
  expect_true(any(grepl('more rows', out)))
  expect_true(any(grepl('N1 = (3, 4), N2 = (2, 5)', out, fixed = TRUE)))
})

test_that("print.bbssr_srr prints a subset without its attributes as a data frame", {
  rr <- StratRR(c(3, 4), c(2, 5), 0.05, 'Strat-Fisher', 'two.sided', 'central')
  # Selecting columns removes the attributes; selecting rows keeps them
  sub <- rr[, names(rr)]
  expect_null(attr(sub, 'Test'))
  out <- capture.output(print(sub))
  expect_false(any(grepl('Rejection region of a stratified test', out)))
  expect_true(any(grepl('S.min', out, fixed = TRUE)))
})

test_that("print.bbssr_srr shows the two-sided method of the unstratified Fisher test", {
  out <- capture.output(print(StratRR(c(3, 4), c(2, 5), 0.05, 'Unstrat-Fisher', 'two.sided',
                                      'blaker')))
  expect_true(any(grepl('Test             : Unstrat-Fisher', out, fixed = TRUE)))
  expect_true(any(grepl('Two-sided method : blaker', out, fixed = TRUE)))
  out <- capture.output(print(StratRR(c(3, 4), c(2, 5), 0.05, 'Unstrat-Chisq',
                                      'two.sided')))
  expect_false(any(grepl('Two-sided method', out, fixed = TRUE)))
})
