cells <- cbind(c(5, 6, 3, 0, 2, 4), c(1, 2, 0, 0, 3, 4)) + 1

test_that("PairedRR returns a bbssr_prr object of the right shape", {
  for (tst in c('McNemar', 'McNemar-exact', 'McNemar-midP', 'McNemar-uncond',
                'Paired-Score')) {
    RR <- PairedRR(N = 8, alpha = 0.05, Test = tst, n.grid = 20)
    expect_s3_class(RR, 'bbssr_prr')
    expect_true(is.matrix(RR))
    expect_type(as.vector(RR), 'logical')
    expect_equal(dim(RR), c(9L, 9L))
    expect_equal(dimnames(RR), list(a = as.character(0:8), b = as.character(0:8)))
    expect_equal(attr(RR, 'Test'), tst)
    expect_false(any(RR[outer(0:8, 0:8, '+') > 8]))
    expect_true(all(is.na(attr(RR, 'p.value')[outer(0:8, 0:8, '+') > 8])))
  }
})

test_that("PairedRR agrees with the reference p-values and regions", {
  check <- function(key.p, key.n, ..., tol = 1e-10) {
    RR <- PairedRR(N = 12, alpha = 0.05, ...)
    expect_equal(attr(RR, 'p.value')[cells], key.p, tolerance = tol)
    expect_equal(sum(RR), key.n)
  }
  check(c(0.0512352174298747, 0.0786496035251426, 0.0416322583317752, 0.5,
          0.672639576990712, 0.5), 21, Test = 'McNemar')
  check(c(0.102470434859749, 0.157299207050285, 0.0832645166635504, 1, 0.654720846018577,
          1), 32, Test = 'McNemar', alternative = 'two.sided')
  check(c(0.109375, 0.14453125, 0.125, 1, 0.8125, 0.63671875), 15, Test = 'McNemar-exact')
  check(c(0.21875, 0.2890625, 0.25, 1, 1, 1), 24, Test = 'McNemar-exact',
        alternative = 'two.sided')
  check(c(0.0625, 0.08984375, 0.0625, 0.5, 0.65625, 0.5), 19, Test = 'McNemar-midP')
  check(c(0.125, 0.1796875, 0.125, 0.5, 0.6875, 0.86328125), 30, Test = 'McNemar-midP',
        alternative = 'two.sided')
  check(c(0.072998046875, 0.113135472540419, 0.072998046875, 1, 1, 1), 18,
        Test = 'McNemar-uncond', ref.pvalue = TRUE, tol = 1e-9)
  check(c(0.14599609375, 0.226270945080838, 0.14599609375, 1, 0.7744140625, 1), 30,
        Test = 'McNemar-uncond', alternative = 'two.sided', ref.pvalue = TRUE, tol = 1e-9)
  check(c(0.0682777946325789, 0.114135472540419, 0.0539247380865849, 1, 0.767453288352893,
          0.635274702815146), 18, Test = 'McNemar-uncond', bb.gamma = 0.001,
        ref.pvalue = TRUE, tol = 1e-9)
  check(c(0.0194153301870173, 0.0340185303947936, 0.0153767806296373, 0.124106539494962,
          0.464034985075496, 0.335099118237478), 27, Test = 'Paired-Score', margin = 0.1)
  check(c(0.0788013610581393, 0.112936259332302, 0.0652850090578681, 0.786611631733508,
          0.761143596853401, 0.584075812604425), 18, Test = 'Paired-Score',
        margin = -0.05)
  check(c(0.88243705602679, 0.843059597448514, 0.899587438652427, 0.124106539494962,
          0.166537180624793, 0.335099118237478), 27, Test = 'Paired-Score',
        alternative = 'less', margin = 0.1)
  check(c(0.829587827230053, 0.788380835652461, 0.834012266458632, 0.124106539494962,
          0.464034985075496, 0.335099118237478), 0, Test = 'Paired-Score',
        alternative = 'equivalence', margin = c(0.1, 0.15))
})

test_that("PairedRR reproduces the asymptotic McNemar p-values of Fagerland et al.", {
  # Table 6: two-sided p-values 0.0339 (N = 21, n12 = 1, n21 = 7) and 0.0330 (N = 161,
  # n12 = 6, n21 = 16)
  expect_equal(round(attr(PairedRR(21, 0.05, 'McNemar', 'two.sided'), 'p.value')[2, 8], 4),
               0.0339)
  expect_equal(round(attr(PairedRR(161, 0.05, 'McNemar', 'two.sided'),
                          'p.value')[7, 17], 4), 0.0330)
})

test_that("PairedRR reproduces the example of Berger and Sidik (2003, Table 3)", {
  # x01 = 4, x10 = 11, n = 30: CB 0.0592, M 0.0354 (Z = 1.8074), Z 0.0494, Z_C 0.0382
  p <- function(...) attr(PairedRR(30, 0.05, ...), 'p.value')[12, 5]
  expect_equal(round(p('McNemar-exact'), 4), 0.0592)
  expect_equal(round(p('McNemar'), 4), 0.0354)
  expect_equal(round(p('McNemar-uncond', ref.pvalue = TRUE), 4), 0.0494)
  expect_equal(round(p('McNemar-uncond', ref.pvalue = TRUE, bb.gamma = 5e-4), 4), 0.0382)
  # Three of the four exact unconditional tests reject at 0.05, the conditional test not
  expect_true(PairedRR(30, 0.05, 'McNemar-uncond', ref.pvalue = TRUE)[12, 5])
  expect_false(PairedRR(30, 0.05, 'McNemar-exact')[12, 5])
})

test_that("PairedRR rejects only in the direction of the alternative", {
  d <- outer(0:10, 0:10, '-')
  for (tst in c('McNemar', 'McNemar-exact', 'McNemar-midP', 'McNemar-uncond')) {
    expect_true(all(d[PairedRR(10, 0.1, tst, n.grid = 30)] > 0), info = tst)
    expect_true(all(d[PairedRR(10, 0.1, tst, 'less', n.grid = 30)] < 0), info = tst)
  }
})

test_that("PairedRR stores the settings", {
  RR <- PairedRR(15, 0.025, 'Paired-Score', 'equivalence', margin = 0.2)
  expect_equal(attr(RR, 'margin'), c(0.2, 0.2))
  expect_equal(attr(RR, 'N'), 15L)
  expect_equal(attr(RR, 'alpha'), 0.025)
  expect_equal(attr(RR, 'alternative'), 'equivalence')
  expect_false(attr(RR, 'ref.pvalue'))
})
