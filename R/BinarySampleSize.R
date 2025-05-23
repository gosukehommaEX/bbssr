########################################################## README ##########################################################
## The code was written by Gosuke Homma (my.name.is.gosuke@gmail.com).
#
## The "BinarySampleSize" function aims to calculate the required sample size for two-arm trials with binary endpoints.
#  User can specify the following five tests:
#  (1) The one-sided Pearson chi-squared test (Chisq)
#  (2) The Fisher exact test (Fisher)
#  (3) The Fisher mid-p test (Fisher-midP)
#  (4) The Z-pooled exact unconditional test (Z-pool)
#  (5) The Boschloo exact unconditional test (Boschloo)
#  Note: The function only covers one-sided tests.
#
## BinarySampleSize has the following arguments.
#  p1:         true probability of responders for group 1
#  p2:         true probability of responders for group 2
#  r:          allocation ratio to group 1 (i.e., an allocation ratio of group 1:group 2 = r:1, r > 0)
#  alpha:      one-sided level of significance
#  tar.power:  target power
#  Test:       type of statistical tests (user can specify 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool' or 'Boschloo')
## BinarySampleSize returns the R data.frame format including the following values. 
#  p1:         true probability of responders for group 1
#  p2:         true probability of responders for group 2
#  r:          allocation ratio to group 1
#  alpha:      one-sided level of significance
#  tar.power:  target power
#  Test:       name of the statistical test
#  Power:      calculated power
#  N1:         required sample size of group 1
#  N2:         required sample size of group 2
#  N:          total required sample size
######################################################## How to use ########################################################
# library(Exact)
# library(exact2x2)
# # (Note) the bottom results are powers obtained from other R packages for validation purposes.
# (1) The one-sided Pearson chi-squared test (Chisq)
# BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 2, alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
#    p1  p2 r alpha tar.power  Test     Power  N1 N2   N
# 1 0.4 0.2 2 0.025       0.8 Chisq 0.8028152 124 62 186
# power.exact.test(0.4, 0.2, 124, 62, method = 'pearson chisq', 'greater', 0.025)$power
# [1] 0.8028152
# power.exact.test(0.4, 0.2, 122, 61, method = 'pearson chisq', 'greater', 0.025)$power
# [1] 0.7938549
#
# (2) The Fisher exact test (Fisher)
# BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 3, alpha = 0.025, tar.power = 0.9, Test = 'Fisher')
#    p1  p2 r alpha tar.power   Test     Power  N1 N2   N
# 1 0.5 0.2 3 0.025       0.9 Fisher 0.9102183 117 39 156
# power.exact.test(0.5, 0.2, 117, 39, method = 'fisher', 'greater', 0.025)$power
# [1] 0.9102183
# power.exact.test(0.5, 0.2, 114, 38, method = 'fisher', 'greater', 0.025)$power
# [1] 0.8996206
#
# (3) The Fisher mid-p test (Fisher-midP)
# BinarySampleSize(p1 = 0.6, p2 = 0.3, r = 2, alpha = 0.025, tar.power = 0.9, Test = 'Fisher-midP')
#    p1  p2 r alpha tar.power        Test     Power N1 N2   N
# 1 0.6 0.3 2 0.025       0.9 Fisher-midP 0.9020784 84 42 126
# Power2x2(84, 42, 0.6, 0.3, alpha = 0.025, pvalFunc=
#            function(x1, n1, x2, n2){
#              fisher.exact(matrix(c(x1, n1 - x1, x2, n2 - x2), 2), alt = 'greater', midp = TRUE)$p.value
#            }
# )
# [1] 0.9020784
# Power2x2(82, 41, 0.6, 0.3, alpha = 0.025, pvalFunc=
#            function(x1, n1, x2, n2){
#              fisher.exact(matrix(c(x1, n1 - x1, x2, n2 - x2), 2), alt = 'greater', midp = TRUE)$p.value
#            }
# )
# [1] 0.8906476
#
# (4) The Z-pooled exact unconditional test (Z-pool)
# BinarySampleSize(p1 = 0.3, p2 = 0.2, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Z-pool')
#    p1  p2 r alpha tar.power   Test     Power  N1  N2   N
# 1 0.3 0.2 1 0.025       0.8 Z-pool 0.8007911 297 297 594
# power.exact.test(0.3, 0.2, 297, 297, method = 'z-pooled', 'greater', 0.025)$power
# [1] 0.8007911
# power.exact.test(0.3, 0.2, 296, 296, method = 'z-pooled', 'greater', 0.025)$power
# [1] 0.799167
#
# (5) The Boschloo exact unconditional test (Boschloo)
# BinarySampleSize(p1 = 0.7, p2 = 0.2, r = 4, alpha = 0.025, tar.power = 0.8, Test = 'Boschloo')
#    p1  p2 r alpha tar.power     Test     Power N1 N2  N
# 1 0.7 0.2 4 0.025       0.8 Boschloo 0.8374957 40 10 50
# power.exact.test(0.7, 0.2, 40, 10, method = 'boschloo', 'greater', 0.025)$power
# [1] 0.8374957
# power.exact.test(0.7, 0.2, 36, 9, method = 'boschloo', 'greater', 0.025)$power
# [1] 0.7998365
############################################################################################################################
source('BinaryPower.R')
BinarySampleSize = function(p1, p2, r, alpha, tar.power, Test) {
  # Step 0 (calculate the required sample size for the one-sided Pearson chi-squared test based on the normal approximation)
  p = (r * p1 + p2) / (1 + r)
  init_N2 = '*'(
    (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha) * sqrt(p * (1 - p)) + qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  # Step 1 (power calculation given initial sample size)
  N2 = ceiling(init_N2)
  N1 = ceiling(r * N2)
  Power = BinaryPower(p1, p2, N1, N2, alpha, Test)
  # Step 2 (sample size calculation via a grid search algorithm)
  if(Power %>=% tar.power) {
    while(Power %>=% tar.power) {
      N2 = N2 - 1
      N1 = ceiling(r * N2)
      Power = BinaryPower(p1, p2, N1, N2, alpha, Test)
    }
    N2 = N2 + 1
  } else {
    while(Power %<<% tar.power) {
      N2 = N2 + 1
      N1 = ceiling(r * N2)
      Power = BinaryPower(p1, p2, N1, N2, alpha, Test)
    }
  }
  # Step 3 (determine the final sample size)
  N1 = ceiling(r * N2)
  N = N1 + N2
  Power = BinaryPower(p1, p2, N1, N2, alpha, Test)
  # Return result
  result = data.frame(p1, p2, r, alpha, tar.power, Test, Power, N1, N2, N)
  return(result)
}