tabs <- list(x1 = rbind(c(4, 7), c(2, 3), c(6, 9), c(3, 5), c(2, 9)),
             x2 = rbind(c(1, 2), c(5, 4), c(0, 5), c(3, 2), c(1, 2)))
pv <- function(Test, alternative, tsmethod = 'minlike') {
  a <- check_strat_args(c(6, 9), c(7, 5), 0.05, Test, alternative, tsmethod)
  n <- nrow(tabs$x1)
  strat_pvalue_cases(matrix(a$N1, n, 2, byrow = TRUE), matrix(a$N2, n, 2, byrow = TRUE),
                     matrix(as.integer(tabs$x1), n), matrix(as.integer(tabs$x2), n),
                     a$type, a$alt, a$ts)
}

test_that("strat_pvalue_cases agrees with the reference values of the asymptotic tests", {
  expect_equal(pv('CMH', 'greater'),
               c(0.0112617092803826, 0.980916090262173, 0.000266002752569626,
                 0.289529653097993, 0.010215560681726),
               tolerance = 1e-10)
  expect_equal(pv('CMH', 'less'),
               c(0.988738290719617, 0.019083909737827, 0.99973399724743,
                 0.710470346902007, 0.989784439318274),
               tolerance = 1e-10)
  expect_equal(pv('CMH', 'two.sided'),
               c(0.0225234185607653, 0.0381678194756539, 0.000532005505139252,
                 0.579059306195985, 0.0204311213634519),
               tolerance = 1e-10)
  expect_equal(pv('CMH-cc', 'greater'),
               c(0.0294481502545263, 0.954817528549904, 0.00170580151718186,
                 0.430466290997207, 0.0312536838324566),
               tolerance = 1e-10)
  expect_equal(pv('CMH-cc', 'less'),
               c(0.970551849745474, 0.0451824714500955, 0.998294198482818,
                 0.569533709002793, 0.968746316167543),
               tolerance = 1e-10)
  expect_equal(pv('CMH-cc', 'two.sided'),
               c(0.0588963005090525, 0.0903649429001911, 0.00341160303436371,
                 0.860932581994413, 0.0625073676649133),
               tolerance = 1e-10)
  expect_equal(pv('Cochran', 'greater'),
               c(0.00886695039356853, 0.984399554405416, 0.000155745488383693,
                 0.282119495929975, 0.00798651767764384),
               tolerance = 1e-10)
  expect_equal(pv('Cochran', 'less'),
               c(0.991133049606431, 0.0156004455945835, 0.999844254511616,
                 0.717880504070025, 0.992013482322356),
               tolerance = 1e-10)
  expect_equal(pv('Cochran', 'two.sided'),
               c(0.0177339007871371, 0.031200891189167, 0.000311490976767385,
                 0.56423899185995, 0.0159730353552877),
               tolerance = 1e-10)
  expect_equal(pv('Cochran-cc', 'greater'),
               c(0.0248091477906745, 0.960797370983855, 0.00115357841581769,
                 0.427763691643765, 0.0264368182089764),
               tolerance = 1e-10)
  expect_equal(pv('Cochran-cc', 'less'),
               c(0.975190852209325, 0.0392026290161454, 0.998846421584182,
                 0.572236308356235, 0.973563181791024),
               tolerance = 1e-10)
  expect_equal(pv('Cochran-cc', 'two.sided'),
               c(0.0496182955813491, 0.0784052580322908, 0.00230715683163538,
                 0.85552738328753, 0.0528736364179528),
               tolerance = 1e-10)
})

# The last table separates the two-sided conventions minlike and blaker
test_that("strat_pvalue_cases agrees with the reference values of the exact test", {
  expect_equal(pv('Strat-Fisher', 'greater'),
               c(0.0271488418341565, 0.994531109915725, 0.000582750582750583,
                 0.429360848941269, 0.0292976254514716),
               tolerance = 1e-10)
  expect_equal(pv('Strat-Fisher', 'less'),
               c(0.997016969044941, 0.0414119679853946, 1,
                 0.829007530930608, 0.99807884423269),
               tolerance = 1e-10)
  expect_equal(pv('Strat-Fisher', 'two.sided', 'minlike'),
               c(0.0435799631603827, 0.056416776696497, 0.000582750582750583,
                 0.706165338158345, 0.0292976254514716),
               tolerance = 1e-10)
  expect_equal(pv('Strat-Fisher', 'two.sided', 'central'),
               c(0.054297683668313, 0.0828239359707891, 0.00116550116550117,
                 0.858721697882537, 0.0585952509029432),
               tolerance = 1e-10)
  expect_equal(pv('Strat-Fisher', 'two.sided', 'blaker'),
               c(0.0435799631603827, 0.056416776696497, 0.000582750582750583,
                 0.706165338158345, 0.0575386152309229),
               tolerance = 1e-10)
})

test_that("strat_pvalue_cases agrees with stats::mantelhaen.test", {
  N1 <- c(6, 9)
  N2 <- c(7, 5)
  for (i in seq_len(nrow(tabs$x1))) {
    x1 <- tabs$x1[i, ]
    x2 <- tabs$x2[i, ]
    tab <- array(0, c(2, 2, 2))
    for (j in 1:2) tab[, , j] <- matrix(c(x1[j], x2[j], N1[j] - x1[j], N2[j] - x2[j]), 2)
    for (alt in c('greater', 'less', 'two.sided')) {
      for (cc in c(FALSE, TRUE)) {
        ref <- stats::mantelhaen.test(tab, alternative = alt, correct = cc)$p.value
        a <- check_strat_args(N1, N2, 0.05, if (cc) 'CMH-cc' else 'CMH', alt, 'minlike')
        got <- strat_pvalue_cases(matrix(as.integer(N1), 1), matrix(as.integer(N2), 1),
                                  matrix(as.integer(x1), 1), matrix(as.integer(x2), 1),
                                  a$type, a$alt, a$ts)
        expect_equal(got, ref, tolerance = 1e-12, info = paste(i, alt, cc))
      }
      # The exact test of mantelhaen.test() orders the two-sided outcomes by their
      # probability, as tsmethod = 'minlike'
      ref <- stats::mantelhaen.test(tab, alternative = alt, exact = TRUE)$p.value
      a <- check_strat_args(N1, N2, 0.05, 'Strat-Fisher', alt, 'minlike')
      got <- strat_pvalue_cases(matrix(as.integer(N1), 1), matrix(as.integer(N2), 1),
                                matrix(as.integer(x1), 1), matrix(as.integer(x2), 1),
                                a$type, a$alt, a$ts)
      expect_equal(got, ref, tolerance = 1e-12, info = paste(i, alt, 'exact'))
    }
  }
})

test_that("strat_pvalue_cases reproduces the example of Jung (2014, Table 2)", {
  # Thymosin (group 1) against placebo in three strata: one-sided p-values 0.1563 for the
  # stratified exact test and 0.0760 for the Mantel-Haenszel test without correction
  N1 <- matrix(c(11L, 9L, 8L), 1)
  N2 <- matrix(c(13L, 12L, 10L), 1)
  x1 <- matrix(c(10L, 9L, 8L), 1)
  x2 <- matrix(c(12L, 11L, 7L), 1)
  expect_equal(round(strat_pvalue_cases(N1, N2, x1, x2, 4L, 0L, 0L), 4), 0.1563)
  expect_equal(round(strat_pvalue_cases(N1, N2, x1, x2, 0L, 0L, 0L), 4), 0.0760)
})

test_that("strat_pvalue_cases gives 1 when the variance is 0 and checks its input", {
  one <- function(x1, x2, type) {
    strat_pvalue_cases(matrix(c(3L, 4L), 1), matrix(c(2L, 5L), 1), matrix(x1, 1),
                       matrix(x2, 1), type, 0L, 0L)
  }
  expect_equal(one(c(0L, 0L), c(0L, 0L), 0L), 1)
  expect_equal(one(c(3L, 4L), c(2L, 5L), 2L), 1)
  expect_error(one(c(4L, 0L), c(0L, 0L), 0L), 'between 0 and the group sizes')
  expect_error(strat_pvalue_cases(matrix(1L, 1, 2), matrix(1L, 1, 3), matrix(0L, 1, 2),
                                  matrix(0L, 1, 2), 0L, 0L, 0L), 'same dimensions')
})

test_that("strat_pvalue_cases agrees with the reference values of the unstratified tests", {
  expect_equal(pv('Unstrat-Chisq', 'greater'),
               c(0.00625122055412587, 0.984345952848196, 0.00029417189438257,
                 0.273295021227503, 0.00625122055412587),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Chisq', 'less'),
               c(0.993748779445874, 0.0156540471518043, 0.999705828105617,
                 0.726704978772497, 0.993748779445874),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Chisq', 'two.sided'),
               c(0.0125024411082517, 0.0313080943036085, 0.000588343788765141,
                 0.546590042455007, 0.0125024411082517),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Chisq-cc', 'greater'),
               c(0.0174260773092047, 0.961266029316397, 0.00137204182392903,
                 0.414761464413403, 0.0174260773092047),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Chisq-cc', 'less'),
               c(0.982573922690795, 0.0387339706836033, 0.998627958176071,
                 0.585238535586597, 0.982573922690795),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Chisq-cc', 'two.sided'),
               c(0.0348521546184094, 0.0774679413672065, 0.00274408364785806,
                 0.829522928826806, 0.0348521546184094),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Fisher', 'greater'),
               c(0.0165320590478754, 0.995231151194269, 0.000891861761426979,
                 0.415304886256562, 0.0165320590478754),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Fisher', 'less'),
               c(0.998439299442126, 0.0377058374837349, 1, 0.838780454973751,
                 0.998439299442126),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Fisher', 'two.sided', 'minlike'),
               c(0.0213009078536067, 0.0542378965316104, 0.000891861761426979,
                 0.703567101897967, 0.0213009078536067),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Fisher', 'two.sided', 'central'),
               c(0.0330641180957509, 0.0754116749674698, 0.00178372352285396,
                 0.830609772513124, 0.0330641180957509),
               tolerance = 1e-10)
  expect_equal(pv('Unstrat-Fisher', 'two.sided', 'blaker'),
               c(0.0213009078536067, 0.0542378965316104, 0.000891861761426979,
                 0.703567101897967, 0.0213009078536067),
               tolerance = 1e-10)
})

test_that("the unstratified tests give the p-values of the table of totals", {
  # Totals of the tables: 15 and 12 patients in the groups
  idx <- cbind(rowSums(tabs$x1) + 1, rowSums(tabs$x2) + 1)
  expect_equal(pv('Unstrat-Chisq', 'greater'),
               stats::pnorm(zstat(15, 12)[idx], lower.tail = FALSE))
  expect_equal(pv('Unstrat-Chisq', 'two.sided'),
               pmin(1, 2 * stats::pnorm(abs(zstat(15, 12)[idx]), lower.tail = FALSE)))
  expect_equal(pv('Unstrat-Fisher', 'greater'),
               fisher_pvalue(15, 12, 'greater', 'minlike', FALSE)[idx])
  for (ts in c('minlike', 'central', 'blaker')) {
    expect_equal(pv('Unstrat-Fisher', 'two.sided', ts),
                 fisher_pvalue(15, 12, 'two.sided', ts, FALSE)[idx])
  }
  # A table whose stratum totals differ but whose totals agree has the same p-value
  a <- check_strat_args(c(6, 9), c(7, 5), 0.05, 'Unstrat-Fisher', 'greater', 'minlike')
  p <- strat_pvalue_cases(rbind(a$N1, a$N1), rbind(a$N2, a$N2), rbind(c(4L, 7L), c(6L, 5L)),
                          rbind(c(1L, 2L), c(0L, 3L)), a$type, a$alt, a$ts)
  expect_identical(p[1], p[2])
})
