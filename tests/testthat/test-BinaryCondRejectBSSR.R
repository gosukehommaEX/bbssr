# Reference values are computed with an independent Python implementation
# (tools/reference/reference_values.py) that sums exact hypergeometric probabilities. The
# first design has Delta.A = 0.3, N1 = N2 = 12, interim sizes 6 and 6, alpha = 0.025,
# target power 0.8, Fisher's exact test and re-estimation by the normal approximation
crp_args <- list(Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
                 alpha = 0.025, tar.power = 0.8, Test = 'Fisher', ss.method = 'standard')

test_that("BinaryCondRejectBSSR agrees with an independent implementation", {
  res <- do.call(BinaryCondRejectBSSR, crp_args)
  expect_s3_class(res, 'bbssr_crp')
  expect_named(res, c('s', 's2', 'N1', 'N2', 'CRP', 'CRP.total'))
  expect_equal(nrow(res), 567)
  # Rows are ordered by s and then by s2, and s2 runs over the second-stage total
  expect_identical(order(res$s, res$s2), seq_len(nrow(res)))
  map <- attr(res, 'reestimation')
  expect_equal(as.vector(table(res$s)), map$N1 + map$N2 - 12 + 1)
  expect_equal(c(sum(res$CRP), sum(res$CRP.total)),
               c(5.8675989615046, 7.48650061067144), tolerance = 1e-10)
  # Outcomes at which the conditional rejection probability exceeds the nominal level
  above <- res[res$CRP > 0.025, ]
  expect_equal(above$s, c(2, 2, 3, 4, 8, 9, 10, 10))
  expect_equal(above$s2, c(3, 6, 15, 21, 43, 37, 30, 33))
  expect_equal(above$CRP,
               c(0.025974025974026, 0.0253599040255932, 0.025881201141577,
                 0.0257019420646433, 0.0257019420646433, 0.025881201141577,
                 0.0253599040255932, 0.025974025974026), tolerance = 1e-10)
  # Fisher's exact test keeps the conditional size given the total at the nominal level
  expect_equal(max(res$CRP.total), 0.0249961144897897, tolerance = 1e-10)
  expect_equal(attr(res, 'exceed'), c(CRP = 8L, CRP.total = 0L))
  m <- attr(res, 'max')
  expect_equal(m$Quantity, c('CRP', 'CRP.total'))
  expect_equal(m$value, c(0.025974025974026, 0.0249961144897897), tolerance = 1e-10)
  # The first of the two outcomes with the largest CRP is reported
  expect_equal(c(m$s[1], m$s2[1]), c(2, 3))
})

test_that("BinaryCondRejectBSSR agrees for the Boschloo test with r = 2", {
  # Same settings with r = 2, N1 = 20, N2 = 10, interim sizes 10 and 5 and the Boschloo
  # test, values from tools/reference/reference_values.py
  res <- BinaryCondRejectBSSR(Delta.A = 0.3, N1 = 20, N2 = 10, n.interim = c(10, 5),
                              r = 2, alpha = 0.025, tar.power = 0.8, Test = 'Boschloo',
                              ss.method = 'standard')
  m <- attr(res, 'max')
  expect_equal(c(nrow(res), m$s[1], m$s2[1], m$s[2], m$s2[2]), c(769, 2, 5, 1, 9))
  expect_equal(c(sum(res$CRP), sum(res$CRP.total)),
               c(12.4601516983984, 16.136283680292), tolerance = 1e-10)
  expect_equal(m$value, c(0.044042913608131, 0.0507470722842913), tolerance = 1e-10)
  # The Boschloo test does not keep the conditional size given the total at the level
  expect_equal(unname(attr(res, 'exceed')), c(173, 265))
})

test_that("the type I error rate from CRP equals that of BinaryTypeIErrorBSSR", {
  theta <- c(0.1, 0.35, 0.6, 0.9)
  cases <- list(
    list(Test = 'Chisq', alternative = 'greater', r = 1, N1 = 15, N2 = 15,
         n.interim = c(7, 7), Delta.A = 0.3),
    list(Test = 'Fisher', alternative = 'less', r = 2, N1 = 20, N2 = 10,
         n.interim = c(8, 4), Delta.A = -0.3),
    list(Test = 'Z-pool', alternative = 'two.sided', r = 1, N1 = 12, N2 = 12,
         n.interim = c(6, 6), Delta.A = 0.3),
    list(Test = 'Boschloo', alternative = 'greater', r = 1, N1 = 12, N2 = 12,
         n.interim = c(6, 5), Delta.A = 0.3)
  )
  for (cs in cases) {
    args <- c(cs, list(alpha = 0.025, tar.power = 0.8, ss.method = 'standard'))
    crp <- do.call(BinaryCondRejectBSSR, c(args, list(theta = theta)))
    tie <- do.call(BinaryTypeIErrorBSSR, c(args, list(theta = theta, refine = FALSE)))
    expect_equal(attr(crp, 'TIE')$TIE, tie$TIE.BSSR, tolerance = 1e-12)
    expect_identical(attr(crp, 'reestimation'), attr(tie, 'reestimation'))
  }
})

test_that("the decomposition by s adds up to the type I error rate", {
  res <- do.call(BinaryCondRejectBSSR, c(crp_args, list(theta = c(0.5, 0.2, 0.2))))
  by.s <- attr(res, 'by.s')
  tie <- attr(res, 'TIE')
  expect_equal(tie$theta, c(0.2, 0.5))
  expect_named(by.s, c('theta', 's', 'N1', 'N2', 'prob.s', 'TIE.s', 'contribution'))
  expect_equal(nrow(by.s), 2 * 13)
  expect_equal(unname(tapply(by.s$prob.s, by.s$theta, sum)), c(1, 1), tolerance = 1e-14)
  expect_equal(unname(tapply(by.s$contribution, by.s$theta, sum)), tie$TIE,
               tolerance = 1e-14)
  expect_true(all(by.s$TIE.s >= 0 & by.s$TIE.s <= 1))
  expect_equal(tie$TIE, c(0.0108006327232362, 0.0158076801034762), tolerance = 1e-10)
})

test_that("without re-estimation CRP averages to CRP.total over the split of the total", {
  res <- do.call(BinaryCondRejectBSSR, c(crp_args, list(N.min = 24, N.max = 24)))
  expect_true(all(res$N1 == 12 & res$N2 == 12))
  expect_gt(nrow(res), 0)
  total <- res$s + res$s2
  # Distribution of the interim count s given the total, among 12 interim and 12
  # second-stage patients
  avg <- rowsum(stats::dhyper(res$s, 12, 12, total) * res$CRP, total, reorder = TRUE)[, 1]
  expect_equal(unname(avg[as.character(total)]), res$CRP.total, tolerance = 1e-12)
  # The two probabilities differ outcome by outcome
  expect_equal(max(abs(res$CRP - res$CRP.total)), 0.019562850663941, tolerance = 1e-8)
})

test_that("ref.pvalue reaches every rejection region and the margin is 0", {
  seen <- ref_pvalue_calls(res <- BinaryCondRejectBSSR(
    Delta.A = 0.5, N1 = 8, N2 = 8, n.interim = c(4, 4), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Z-pool', ss.method = 'standard', ref.pvalue = TRUE
  ))
  # One rejection region for each distinct final sample size
  map <- attr(res, 'reestimation')
  expect_gt(nrow(seen), 1)
  expect_equal(nrow(seen), nrow(unique(map[, c('N1', 'N2')])))
  expect_true(all(seen$ref.pvalue))
  expect_true(all(seen$margin == 0))
  expect_true(attr(res, 'ref.pvalue'))
})

test_that("BinaryCondRejectBSSR validates theta and takes no margin", {
  expect_error(do.call(BinaryCondRejectBSSR, c(crp_args, list(theta = c(0.2, 1.2)))),
               'theta')
  expect_error(do.call(BinaryCondRejectBSSR, c(crp_args, list(margin = 0.1))),
               'unused argument')
})
