pts <- list(p10 = c(0.3, 0.25, 0.05, 0.6), p01 = c(0.1, 0.25, 0, 0.4))

test_that("PairedPower agrees with the reference values", {
  pw <- function(...) PairedPower(pts$p10, pts$p01, N = 25, alpha = 0.05, ...)$Power
  expect_equal(pw(Test = 'McNemar'),
               c(0.501762078611304, 0.0510007386149977, 0.127106495660931,
                 0.273531450144572), tolerance = 1e-10)
  expect_equal(pw(Test = 'McNemar-exact', alternative = 'two.sided'),
               c(0.241948002676632, 0.0247683375474778, 0.00121296132209725,
                 0.154757175258514), tolerance = 1e-10)
  expect_equal(pw(Test = 'McNemar-midP'),
               c(0.455567858074236, 0.0467111657893424, 0.0340906014809905,
                 0.273531450144572), tolerance = 1e-10)
  expect_equal(pw(Test = 'McNemar-cc'),
               c(0.348235072572542, 0.0276578546401005, 0.00716494790258658,
                 0.153551734755933), tolerance = 1e-10)
  expect_equal(pw(Test = 'McNemar-cc', alternative = 'two.sided'),
               c(0.241258649191447, 0.0235977503363145, 0.00121296132209725,
                 0.154757175258514), tolerance = 1e-10)
  expect_equal(pw(Test = 'McNemar-uncond', ref.pvalue = TRUE),
               c(0.425551746919397, 0.0368039092574043, 0.0340906014809905,
                 0.153551734755933), tolerance = 1e-9)
  expect_equal(pw(Test = 'Paired-Score', margin = 0.1),
               c(0.777878646102371, 0.165658340099746, 0.999999999999996,
                 0.424617017671006), tolerance = 1e-10)
  expect_equal(PairedPower(c(0.15, 0.2, 0.05), c(0.15, 0.1, 0), N = 40, alpha = 0.05,
                           Test = 'Paired-Score', alternative = 'equivalence',
                           margin = 0.2)$Power,
               c(0.476798189974026, 0.281649326426356, 0.861850224498939),
               tolerance = 1e-10)
})

test_that("PairedPower reproduces the exact sizes of Berger and Sidik (2003, Table 4)", {
  # One-sided level 0.05; the size is the largest rejection probability on p10 = p01
  size <- function(n, ...) {
    rr <- get_paired_rr(n, 0.05, ..., 'greater', 100, 0, TRUE, 0)
    paired_tie_certify(list(rr), n, list(rep(1, n + 1)), c(0, 1), 'greater', 0)$y
  }
  n <- c(10, 15, 20, 25)
  expect_equal(round(vapply(n, size, numeric(1), 'McNemar-uncond'), 4),
               c(0.0265, 0.0369, 0.0393, 0.0382))
  expect_equal(round(vapply(n, size, numeric(1), 'McNemar'), 4),
               c(0.0652, 0.0592, 0.0577, 0.0539))
  expect_equal(round(vapply(n, size, numeric(1), 'McNemar-exact'), 4),
               c(0.0208, 0.0304, 0.0339, 0.0342))
})

test_that("PairedPower reproduces the average powers of Berger and Sidik (2003, Table 5)", {
  # 100 alternatives p01 = 0.025, ..., 0.475 and p10 = p01 + 0.05, ..., 1 - p01
  alt <- do.call(rbind, lapply(seq(0.025, 0.475, by = 0.05), function(q) {
    data.frame(p01 = q, p10 = seq(q + 0.05, 1 - q + 1e-9, by = 0.05))
  }))
  expect_equal(nrow(alt), 100)
  avg <- function(n, Test, ...) {
    mean(PairedPower(alt$p10, alt$p01, N = n, alpha = 0.05, Test = Test, ...)$Power)
  }
  expect_equal(round(c(avg(10, 'McNemar-exact'), avg(25, 'McNemar-exact')), 3),
               c(0.287, 0.555))
  expect_equal(round(c(avg(10, 'McNemar'), avg(25, 'McNemar')), 3), c(0.439, 0.618))
  expect_equal(round(c(avg(10, 'McNemar-uncond', ref.pvalue = TRUE),
                       avg(25, 'McNemar-uncond', ref.pvalue = TRUE)), 3), c(0.316, 0.577))
})

test_that("PairedPower recycles a single probability", {
  pw <- PairedPower(c(0.3, 0.35), 0.1, N = 20, alpha = 0.05, Test = 'McNemar')
  expect_equal(pw$p01, c(0.1, 0.1))
  expect_equal(pw$Power, PairedPower(c(0.3, 0.35), c(0.1, 0.1), N = 20, alpha = 0.05,
                                     Test = 'McNemar')$Power)
  expect_error(PairedPower(c(0.3, 0.35), c(0.1, 0.1, 0.1), N = 20, alpha = 0.05,
                           Test = 'McNemar'), 'one of them a single value')
})

test_that("PairedPower checks the probabilities and returns the settings", {
  expect_error(PairedPower(0.6, 0.5, N = 10, alpha = 0.05, Test = 'McNemar'),
               'must not exceed 1')
  pw <- PairedPower(0.3, 0.1, N = 10, alpha = 0.05, Test = 'McNemar-e')
  expect_s3_class(pw, 'bbssr_ppower')
  expect_equal(pw$Test, 'McNemar-exact')
  expect_equal(pw$N, 10L)
  expect_named(pw, c('p10', 'p01', 'N', 'alpha', 'Test', 'alternative', 'Power'))
})
