######################################################## README ########################################################
## The code was written by Gosuke Homma (my.name.is.gosuke@gmail.com).
#
## The "BinaryPower" function aims to calculate power for two-arm trials with binary endpoints using an exact test.
#  User can specify the following five tests:
#  (1) The one-sided Pearson chi-squared test (Chisq)
#  (2) The Fisher exact test (Fisher)
#  (3) The Fisher mid-p test (Fisher-midP)
#  (4) The Z-pooled exact unconditional test (Z-pool)
#  (5) The Boschloo exact unconditional test (Boschloo)
#  Note: The function only covers one-sided tests.
#
## BinaryPower has the following arguments.
#  p1:         true probability of responders for group 1 (note: user can specify a vector with different values)
#  p2:         true probability of responders for group 2 (note: user can specify a vector with different values)
#  N1:         sample size for group 1
#  N2:         sample size for group 2
#  alpha:      one-sided level of significance
#  Test:       type of statistical tests (user can specify 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool' or 'Boschloo')
#
## BinaryPower returns the power. If the user specifies vectors of p1 and p2, a vector of powers will be returned.
###################################################### How to use ######################################################
# library(Exact)
# library(exact2x2)
# p1 = c(0.5, 0.6, 0.7, 0.8)
# p2 = c(0.2, 0.2, 0.2, 0.2)
# N1 = 10
# N2 = 40
# alpha = 0.025
# # (Note) the bottom results are powers obtained from other R packages for validation purposes.
# (1) The one-sided Pearson chi-squared test (Chisq)
# BinaryPower(p1, p2, N1, N2, alpha, Test = 'Chisq')
# [1] 0.4847117 0.6971209 0.8664207 0.9638721
# sapply(seq(p1), function(i) {
#   power.exact.test(p1[i], p2[i], N1, N2, method = 'pearson chisq', 'greater', alpha)$power
# })
# [1] 0.4847117 0.6971209 0.8664207 0.9638721
#
# (2) The Fisher exact test (Fisher)
# BinaryPower(p1, p2, N1, N2, alpha, Test = 'Fisher')
# [1] 0.3288391 0.5474891 0.7619361 0.9172846
# sapply(seq(p1), function(i) {
#   power.exact.test(p1[i], p2[i], N1, N2, method = 'fisher', 'greater', alpha)$power
# })
# [1] 0.3288391 0.5474891 0.7619361 0.9172846
#
# (3) The Fisher mid-p test (Fisher-midP)
# BinaryPower(p1, p2, N1, N2, alpha, Test = 'Fisher-midP')
# [1] 0.4057605 0.6275103 0.8230326 0.9475556
# sapply(seq(p1), function(i) {
#   fisher.midP = Power2x2(N1, N2, p1[i], p2[i], alpha, pvalFunc=
#                            function(x1, n1, x2, n2){
#                              fisher.exact(matrix(c(x1, n1 - x1, x2, n2 - x2), 2), alt = 'greater', midp = TRUE)$p.val
#                            }
#   )
# })
# [1] 0.4057605 0.6275103 0.8230326 0.9475556
#
# (4) The Z-pooled exact unconditional test (Z-pool)
# BinaryPower(p1, p2, N1, N2, alpha, Test = 'Z-pool')
# [1] 0.4064897 0.6277024 0.8229305 0.9473531
# sapply(seq(p1), function(i) {
#   power.exact.test(p1[i], p2[i], N1, N2, method = 'z-pooled', 'greater', alpha)$power
# })
# [1] 0.4064897 0.6277024 0.8229305 0.9473531
#
# (5) The Boschloo exact unconditional test (Boschloo)
# BinaryPower(p1, p2, N1, N2, alpha, Test = 'Boschloo')
# [1] 0.4277969 0.6544621 0.8445363 0.9570202
# sapply(seq(p1), function(i) {
#   power.exact.test(p1[i], p2[i], N1, N2, method = 'boschloo', 'greater', alpha)$power
# })
# [1] 0.4277969 0.6544621 0.8445363 0.9570202
########################################################################################################################
source('BinaryRR.R')
BinaryPower = function(p1, p2, N1, N2, alpha, Test) {
  # Check that p1 and p2 are the same length.
  if(length(p1) %!=% length(p2)) stop('x1 and x2 should be the same length')
  # Set rejection region
  RR = BinaryRR(N1, N2, alpha, Test)
  # Return power
  sapply(seq(length(p1)), function(i) sum(dbinom(0:N1, N1, p1[i]) * pbinom(rowSums(RR) - 1, N2, p2[i])))
}