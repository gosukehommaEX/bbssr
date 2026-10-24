scrp <- function() {
  StratCondRejectBSSR(Delta.A = 0.35, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                      alpha = 0.05, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
                      N.max = 36, theta = c(0.3, 0.5))
}

test_that("print.bbssr_scrp prints the settings and the largest probabilities", {
  crp <- scrp()
  out <- capture.output(r <- print(crp))
  expect_identical(r, crp)
  expect_true(any(grepl('for a stratified design', out, fixed = TRUE)))
  expect_true(any(grepl('Interim sizes    : (3, 3) and (3, 3)', out, fixed = TRUE)))
  expect_true(any(grepl('Interim outcomes : 49 pairs (s1, s2)', out, fixed = TRUE)))
  expect_true(any(grepl('Type I error rate over 4 points', out, fixed = TRUE)))
})

test_that("print.bbssr_scrp prints a subset without its attributes as a data frame", {
  crp <- scrp()
  # Selecting columns removes the attributes; selecting rows keeps them
  sub <- crp[, names(crp)]
  expect_null(attr(sub, 'max'))
  out <- capture.output(print(sub))
  expect_false(any(grepl('for a stratified design', out, fixed = TRUE)))
  expect_true(any(grepl('CRP.max', out, fixed = TRUE)))
})
