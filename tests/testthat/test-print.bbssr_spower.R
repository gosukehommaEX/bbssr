test_that("print.bbssr_spower prints the settings and the power", {
  res <- StratPower(c(0.6, 0.45), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'CMH')
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl('Exact power of a stratified test', out)))
  expect_true(any(grepl('0.3033', out, fixed = TRUE)))
})

test_that("print.bbssr_spower prints a subset without its attributes as a data frame", {
  res <- StratPower(c(0.6, 0.45), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'CMH')
  # Selecting columns removes the attributes; selecting rows keeps them
  sub <- res[, names(res)]
  expect_null(attr(sub, 'Test'))
  out <- capture.output(print(sub))
  expect_false(any(grepl('Exact power of a stratified test', out)))
  expect_true(any(grepl('Power', out, fixed = TRUE)))
})

test_that("print.bbssr_spower shows the two-sided method of the unstratified Fisher test", {
  res <- StratPower(c(0.6, 0.45), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'Unstrat-Fisher',
                    'two.sided', 'central')
  out <- capture.output(print(res))
  expect_true(any(grepl('Two-sided method : central', out, fixed = TRUE)))
})

test_that("print.bbssr_spower describes a random allocation", {
  res <- StratPower(c(0.5, 0.8), c(0.2, 0.4), alpha = 0.05, Test = 'Strat-Fisher',
                    alloc = 'group', N = 17, frac = c(0.4, 0.6), r = 3)
  out <- capture.output(print(res))
  expect_true(any(grepl('stratum sizes (6, 11), random group sizes', out, fixed = TRUE)))
  expect_true(any(grepl('Total size       : N = 17, fractions (0.4, 0.6), ratio 3', out,
                        fixed = TRUE)))
  expect_false(any(grepl('Group sizes', out, fixed = TRUE)))
  res <- StratPower(c(0.5, 0.8), c(0.2, 0.4), alpha = 0.05, Test = 'CMH', alloc = 'stratum',
                    N = 16, frac = c(0.5, 0.5))
  out <- capture.output(print(res))
  expect_true(any(grepl('permuted blocks of 4 within the strata', out, fixed = TRUE)))
  expect_true(any(grepl('Left out', out, fixed = TRUE)))
})
