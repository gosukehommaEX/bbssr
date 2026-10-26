pbssr <- function(...) {
  args <- list(psi = c(0.3, 0.5), Delta.T = 0.25, Delta.A = 0.25, N = 30, n.interim = 10,
               alpha = 0.05, tar.power = 0.8, Test = 'McNemar', ss.method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  do.call(PairedPowerBSSR, args)
}

test_that("PairedPowerBSSR agrees with the reference values", {
  get <- function(res) c(rbind(res$power.BSSR, res$E.N))
  expect_equal(get(pbssr()), c(0.842667723200926, 31.6893289541, 0.779391115235398,
                               48.3798828125), tolerance = 1e-10)
  expect_equal(get(pbssr(ss.method = 'exact')),
               c(0.759431283045512, 27.2010519196, 0.754755057669941, 46.41796875),
               tolerance = 1e-10)
  expect_equal(get(pbssr(Test = 'McNemar-exact', alternative = 'two.sided',
                         ss.method = 'conditional.power', N.max = 45)),
               c(0.787666291589832, 38.3941968884, 0.580758613628438, 43.35546875),
               tolerance = 1e-10)
  expect_equal(get(pbssr(Test = 'McNemar-cc', ss.method = 'exact')),
               c(0.77340819319093, 35.6456938805, 0.757314659877442, 53.8486328125),
               tolerance = 1e-10)
  expect_equal(get(pbssr(Delta.T = 0, Delta.A = 0, Test = 'Paired-Score', margin = 0.2)),
               c(0.814580182356671, 50.9347553717, 0.795139257151903, 77.9501953125),
               tolerance = 1e-10)
})

test_that("PairedPowerBSSR agrees with enumeration of both stages", {
  # A large effect keeps the final numbers of pairs small for the enumeration
  res <- pbssr(psi = 0.5, Delta.T = 0.4, Delta.A = 0.4, N = 20, n.interim = 6)
  expect_lte(max(attr(res, 'reestimation')$N), 40)
  map <- attr(res, 'reestimation')
  Ns <- sort(unique(map$N))
  regions <- lapply(Ns, function(n) {
    get_paired_rr(n, 0.05, 'McNemar', 'greater', 100, 0, FALSE, 0)
  })
  names(regions) <- Ns
  expect_equal(res$power.BSSR, paired_bssr_ref(map$N, regions, 6, res$p10, res$p01),
               tolerance = 1e-12)
})

test_that("PairedPowerBSSR returns the fixed design and the distribution of N", {
  res <- pbssr(psi = c(0.1, 0.3, 1.2))
  # psi = 0.1 is below Delta.T and psi = 1.2 above 1
  expect_equal(res$psi, 0.3)
  expect_equal(res$power.TRAD, PairedPower(0.275, 0.025, 30, 0.05, 'McNemar')$Power)
  d <- attr(res, 'N.dist')
  expect_equal(sum(d$prob), 1)
  expect_equal(sum(d$N * d$prob), res$E.N)
  expect_s3_class(res, 'bbssr_ppowerbssr')
  expect_equal(attr(res, 'n.interim'), 10L)
  expect_error(pbssr(psi = 0.1), 'between |Delta.T| and 1', fixed = TRUE)
  expect_error(pbssr(Delta.T = 1), 'Delta.T')
})

test_that("the exact conditional McNemar test keeps its level after re-estimation", {
  res <- pbssr(psi = c(0.1, 0.3, 0.6, 0.9), Delta.T = 0, Test = 'McNemar-exact',
               alternative = 'two.sided', ss.method = 'conditional.power', N.max = 45)
  expect_true(all(res$power.BSSR <= 0.05))
})
