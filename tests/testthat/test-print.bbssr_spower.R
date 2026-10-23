test_that("print.bbssr_spower prints the settings and the power", {
  res <- StratPower(c(0.6, 0.45), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'CMH')
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl('Exact power of a stratified test', out)))
  expect_true(any(grepl('0.3033', out, fixed = TRUE)))
})

test_that("print.bbssr_spower prints a subset without its attributes as a data frame", {
  res <- StratPower(c(0.6, 0.45), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'CMH')
  sub <- res[1, ]
  expect_null(attr(sub, 'Test'))
  out <- capture.output(print(sub))
  expect_false(any(grepl('Exact power of a stratified test', out)))
  expect_true(any(grepl('Power', out, fixed = TRUE)))
})
