test_that("print.bbssr_prr reports the settings and the map", {
  RR <- PairedRR(N = 6, alpha = 0.1, Test = 'McNemar-uncond', bb.gamma = 0.001)
  out <- utils::capture.output(ret <- print(RR))
  expect_identical(ret, RR)
  expect_true(any(grepl('Pairs         : N = 6', out, fixed = TRUE)))
  expect_true(any(grepl('Berger-Boos   : 0.001', out, fixed = TRUE)))
  expect_true(any(grepl(sprintf('Rejected outcomes: %d of 28', sum(RR)), out, fixed = TRUE)))
  expect_true(any(grepl('a=6', out, fixed = TRUE)))
})

test_that("print.bbssr_prr suppresses the map for many pairs", {
  RR <- PairedRR(N = 30, alpha = 0.05, Test = 'Paired-Score', margin = 0.1)
  out <- utils::capture.output(print(RR))
  expect_true(any(grepl('Outcome map suppressed', out)))
  expect_true(any(grepl('Margin        : 0.1', out, fixed = TRUE)))
  out <- utils::capture.output(print(RR, show.map = TRUE))
  expect_true(any(grepl('b=30', out, fixed = TRUE)))
})
