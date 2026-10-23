count_rejected <- function(rr) {
  sum(pmax(0, rr$lower - rr$S.min + 1), na.rm = TRUE) +
    sum(pmax(0, rr$S.max - rr$upper + 1), na.rm = TRUE)
}

test_that("StratRR agrees with the reference numbers of rejected outcomes", {
  n <- function(Test, alternative) {
    count_rejected(StratRR(c(5, 7), c(6, 4), 0.1, Test, alternative, 'blaker'))
  }
  expect_equal(c(n('CMH', 'greater'), n('CMH', 'two.sided')), c(253, 396))
  expect_equal(c(n('CMH-cc', 'greater'), n('CMH-cc', 'two.sided')), c(179, 266))
  expect_equal(c(n('Cochran', 'greater'), n('Cochran', 'two.sided')), c(256, 420))
  expect_equal(c(n('Cochran-cc', 'greater'), n('Cochran-cc', 'two.sided')), c(187, 294))
  expect_equal(c(n('Strat-Fisher', 'greater'), n('Strat-Fisher', 'two.sided')),
               c(186, 342))
})

test_that("StratRR lists the stratum totals with the first stratum varying fastest", {
  rr <- StratRR(c(3, 4), c(2, 5), 0.05, 'CMH', 'two.sided')
  expect_s3_class(rr, 'bbssr_srr')
  expect_named(rr, c('z.S1', 'z.S2', 'S.min', 'S.max', 'lower', 'upper'))
  expect_identical(nrow(rr), 6L * 10L)
  expect_identical(rr$z.S1[1:7], c(0:5, 0L))
  expect_identical(rr$z.S2[1:7], c(rep(0L, 6), 1L))
  expect_identical(rr$S.min, as.integer(pmax(0, rr$z.S1 - 2) + pmax(0, rr$z.S2 - 5)))
  expect_identical(rr$S.max, as.integer(pmin(rr$z.S1, 3) + pmin(rr$z.S2, 4)))
  ok <- !is.na(rr$upper)
  expect_true(all(rr$upper[ok] >= rr$S.min[ok] & rr$upper[ok] <= rr$S.max[ok]))
  ok <- !is.na(rr$lower)
  expect_true(all(rr$lower[ok] >= rr$S.min[ok] & rr$lower[ok] <= rr$S.max[ok]))
  expect_identical(attr(rr, 'Test'), 'CMH')
})

test_that("StratRR rejects in one tail for a one-sided alternative", {
  up <- StratRR(c(6, 6), c(6, 6), 0.05, 'Strat-Fisher', 'greater')
  down <- StratRR(c(6, 6), c(6, 6), 0.05, 'Strat-Fisher', 'less')
  expect_true(all(is.na(up$lower)) && any(!is.na(up$upper)))
  expect_true(all(is.na(down$upper)) && any(!is.na(down$lower)))
  # With equal groups the alternative 'less' mirrors 'greater': S becomes S.max + S.min - S
  expect_identical(down$S.max + down$S.min - down$lower, up$upper)
})

test_that("StratRR refuses designs with too many combinations", {
  expect_error(StratRR(rep(10, 6), rep(10, 6), 0.05, 'CMH'), 'more than 1e6')
})
