######################################################## README ########################################################
## The code was written by Gosuke Homma (my.name.is.gosuke@gmail.com).
#
## The "BinaryRR" function aims to provide a rejection region (RR) for two-arm trials with binary endpoints.
#  User can specify the following five tests:
#  (1) The one-sided Pearson chi-squared test (Chisq)
#  (2) The Fisher exact test (Fisher)
#  (3) The Fisher mid-p test (Fisher-midP)
#  (4) The Z-pooled exact unconditional test (Z-pool)
#  (5) The Boschloo exact unconditional test (Boschloo)
#  Note: The function only covers one-sided tests.
#
## BinaryRR has the following arguments.
#  N1:         sample size for group 1
#  N2:         sample size for group 2
#  alpha:      one-sided level of significance
#  Test:       type of the statistical test (user can specify 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool' or 'Boschloo')
#
## BinaryRR returns the RR in R matrix format.
###################################################### How to use ######################################################
# N1 = 20
# N2 = 10
# alpha = 0.025
# Test = 'Boschloo'
# BinaryRR(N1, N2, alpha, Test)
#        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11]
# [1,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [2,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [3,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [4,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [5,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [6,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [7,]  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [8,]   TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [9,]   TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [10,]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [11,]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [12,]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [13,]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [14,]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [15,]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [16,]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [17,]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
# [18,]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
# [19,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
# [20,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
# [21,]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE
########################################################################################################################
library(fpCompare)
BinaryRR = function(N1, N2, alpha, Test) {
  if((Test == 'Chisq') | (Test == 'Z-pool')) {
    # Test statistics for the Chisq over all combinations of i(=0,...,N1) and j(=0,...,N2)
    Z.ij = '/'(
      outer(N2 * (0:N1), N1 * (0:N2), '-') / (N1 * N2),
      sqrt(outer(0:N1, 0:N2, '+') / (N1 * N2) * (1 - outer(0:N1, 0:N2, '+') / (N1 + N2)))
    )
    Z.ij[is.na(Z.ij)] = 0
    if(Test == 'Chisq') {
      # p-values for the Chisq
      p.val = pnorm(Z.ij, lower.tail = FALSE)
    } else if(Test == 'Z-pool') {
      # Since zero and negative values of test statistics must not be statistically significant, they are omitted
      Z.ij.posi = Z.ij[Z.ij %>>% 0]
      order.Z.ij.posi = order(Z.ij.posi, decreasing = TRUE)
      i = (row(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]
      j = (col(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]
      # Calculate P_H0(X1 = i, X2 = j | theta)
      uniq.i = sort(unique(i))
      dbinom.i = sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, N1, theta))
      uniq.j = sort(unique(j))
      dbinom.j = sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, N2, theta))
      P_H0 = dbinom.i[i, ] * dbinom.j[j + 1, ]
      # p-values for all possible values of theta
      p.ij = apply(apply(P_H0, 2, cumsum), 1, max)
      # p-values for the Z-pool
      p.val = 1 ^ Z.ij
      p.val[cbind(i + 1, j + 1)] = p.ij
    }
  } else if(Test == 'Fisher-midP') {
    # p-values for the Fisher-midP
    p.val = outer(0:N1, 0:N2, function(i, j) phyper(i, N1, N2, i + j, lower.tail = FALSE) + 0.5 * dhyper(i, N1, N2, i + j))
  } else if((Test == 'Fisher') | (Test == 'Boschloo')) {
    # p-values for the Fisher over all combinations of i(=0,...,N1) and j(=0,...,N2)
    p.fisher = outer(0:N1, 0:N2, function(i, j) phyper(i - 1, N1, N2, i + j, lower.tail = FALSE))
    if(Test == 'Fisher') {
      # p-values for the Fisher
      p.val = p.fisher
    } else if(Test == 'Boschloo') {
      # Since p-values satisfying hat{p}_{1} - hat{p}_{2} <= 0 must not be statistically significant, they are omitted
      p.max.boschloo = min(p.fisher[(outer(N2 * (0:N1), N1 * (0:N2), '-') / (N1 * N2)) %<=% 0])
      p.fisher.posi = p.fisher[p.fisher %<<% p.max.boschloo]
      order.p.fisher.posi = order(p.fisher.posi, decreasing = FALSE)
      i = (row(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]
      j = (col(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]
      # Calculate P_H0(X1 = i, X2 = j | theta)
      uniq.i = sort(unique(i))
      dbinom.i = sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, N1, theta))
      uniq.j = sort(unique(j))
      dbinom.j = sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, N2, theta))
      P_H0 = dbinom.i[i - min(uniq.i) + 1, ] * dbinom.j[j + 1, ]
      # p-values for all possible values of theta
      p.ij = apply(apply(P_H0, 2, cumsum), 1, max)
      # p-values for the Boschloo
      p.val = 1 ^ p.fisher
      p.val[cbind(i + 1, j + 1)] = p.ij
    }
  }
  # Rejection region
  RR = (p.val %<<% alpha)
  # Return RR
  return(RR)
}