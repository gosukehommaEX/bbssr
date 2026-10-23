test_that("print.bbssr_pbssr reports the re-estimation", {
  res <- PairedBSSR(n1 = 10, D1 = 4, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                    Test = 'McNemar', ss.method = 'standard')
  out <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl('n1 = 10, discordant D1 = 4', out, fixed = TRUE)))
  expect_true(any(grepl('Final N          : 38', out, fixed = TRUE)))
  expect_true(any(grepl('Conditional power', out, fixed = TRUE)))
})

test_that("print.bbssr_pbssr prints a subset without its attributes as a data frame", {
  res <- PairedBSSR(n1 = 10, D1 = 4, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                    Test = 'McNemar', ss.method = 'standard')
  sub <- res[1, ]
  expect_null(attr(sub, 'ss.method'))
  out <- utils::capture.output(print(sub))
  expect_false(any(grepl('paired interim data', out, fixed = TRUE)))
  expect_true(any(grepl('psi.hat', out, fixed = TRUE)))
})
