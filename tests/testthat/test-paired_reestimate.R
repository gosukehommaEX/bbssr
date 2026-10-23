re_args <- list(n1 = 10L, D1 = 0:10, Delta.A = 0.25, N = 30, tar.power = 0.8,
                restricted = FALSE, alternative = 'greater', n.grid = 100, bb.gamma = 0,
                ss.method = 'standard', ss.Test = 'McNemar', ss.alpha = 0.05, N.min = NULL,
                N.max = NULL, ref.pvalue = FALSE, margin = 0, search = 'crossing',
                search.limit = c(2, 50))
re <- function(...) {
  args <- re_args
  new <- list(...)
  args[names(new)] <- new
  do.call(paired_reestimate, args)
}

test_that("paired_reestimate agrees with the reference final numbers of pairs", {
  expect_equal(re()$N, c(23, 23, 23, 28, 38, 48, 58, 68, 78, 87, 97))
  expect_equal(re(ss.method = 'exact')$N, c(16, 16, 16, 24, 35, 47, 57, 68, 77, 87, 97))
  expect_equal(re(ss.method = 'conditional.power', ss.Test = 'McNemar-exact',
                  alternative = 'two.sided', N.max = 45)$N,
               c(41, 36, 31, 37, 45, 45, 45, 45, 45, 45, 45))
  expect_equal(re(Delta.A = 0, ss.Test = 'Paired-Score', margin = 0.2)$N,
               c(11, 26, 38, 50, 64, 77, 92, 106, 121, 136, 151))
})

test_that("paired_reestimate raises the estimate to the assumed difference", {
  r <- re()
  expect_equal(r$psi.hat, (0:10) / 10)
  expect_equal(r$psi.plan, pmax((0:10) / 10, 0.25))
  expect_equal(r$p10.plan - r$p01.plan, rep(0.25, 11))
  expect_equal(r$p10.plan + r$p01.plan, r$psi.plan)
  expect_equal(r$N.re, r$N)
  expect_equal(r$n.raw, paired_raw_n(r$p10.plan, r$p01.plan, 0.05, 0.8, 'greater', 0))
})

test_that("paired_reestimate applies the bounds", {
  r <- re(restricted = TRUE, N.max = 60)
  expect_equal(r$N, pmin(60, pmax(30, re()$N)))
  expect_equal(re(N.min = 40)$N, pmax(40, re()$N))
  expect_error(re(N.min = 50, N.max = 40), 'N.max is smaller')
  expect_error(re(ss.method = 'conditional.power'), 'requires N.max')
  expect_error(re(restricted = TRUE, N = NULL), 'N must be supplied')
  expect_error(re(Delta.A = -0.1), "Delta.A must exceed 0")
  expect_error(re(ss.Test = 'McNemar', margin = 0.1, Delta.A = 0.2), 'Paired-Score')
})

test_that("the conditional power rule takes the smallest number attaining the target", {
  r <- re(ss.method = 'conditional.power', ss.Test = 'McNemar-exact',
          alternative = 'two.sided', N.max = 45)
  cp <- function(N, i) {
    rr <- get_paired_rr(N, 0.05, 'McNemar-exact', 'two.sided', 100, 0, FALSE, 0)
    paired_cp(rr, 10L, r$D1[i], r$psi.plan[i], r$p10.plan[i] / r$psi.plan[i])
  }
  for (i in c(1, 3, 4)) {
    expect_gte(cp(r$N[i], i), 0.8)
    expect_true(all(vapply(10:(r$N[i] - 1), cp, numeric(1), i = i) < 0.8))
  }
  # The upper bound when the target is out of reach
  expect_lt(cp(45, 11), 0.8)
})
