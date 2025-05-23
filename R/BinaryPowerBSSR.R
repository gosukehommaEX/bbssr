######################################################## README ########################################################
## The code was written by Gosuke Homma (my.name.is.gosuke@gmail.com).
#
## The "BinaryPowerBSSR" function aims to calculate the power for two-arm trials with binary endpoints
#  when blinded sample size re-estimation (BSSR) is implemented.
#  User can specify the following five tests:
#  (1) The one-sided Pearson chi-squared test (Chisq)
#  (2) The Fisher exact test (Fisher)
#  (3) The Fisher mid-p test (Fisher-midP)
#  (4) The Z-pooled exact unconditional test (Z-pool)
#  (5) The Boschloo exact unconditional test (Boschloo)
#  Note: The function only covers one-sided tests.
#
## BinaryPowerBSSR has the following arguments.
#  asmd.p1:    an assumed proportion of responders for group 1
#  asmd.p2:    an assumed proportion of responders for group 2
#  p:          pooled proportion of responders from both groups (user can specify multiple values as the vector)
#  Delta.A:    an assumed treatment effect of a risk difference
#  Delta.T:    a true treatment effect of a risk difference
#  N1:         initial sample size of group 1
#  N2:         initial sample size of group 2
#  omega:      fraction of sample size used for an interim analysis (i.e., for BSSR)
#  r:          allocation ratio to group 1
#  alpha:      one-sided level of significance
#  tar.power:  target power
#  Test:       type of statistical tests (user can specify 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool' or 'Boschloo')
#  restricted: specify restricted design or unrestricted design (if restricted = TRUE, the restricted design is chosen)
#  weighted:   specify weighted approach (if weighted = TRUE, the weighted approach is chosen)
#
## BinaryPowerBSSR returns the R data.frame format including the following values.
#  p1:         true probability of responders for group 1
#  p2:         true probability of responders for group 2
#  p:          true probability of pooled responders from both groups
#  power.BSSR: power for BSSR
#  power.TRAD: powers for traditional design
###################################################### How to use ######################################################
# BinaryPowerBSSR(
#   asmd.p1 = 0.45,
#   asmd.p2 = 0.09,
#   p = seq(0, 1, by = 0.01),
#   Delta.A = 0.36,
#   Delta.T = 0.36,
#   N1 = 24,
#   N2 = 24,
#   omega = 0.5,
#   r = 1,
#   alpha = 0.025,
#   tar.power = 0.8,
#   Test = 'Z-pool',
#   restricted = FALSE,
#   weighted = TRUE
# )
# 
#           p1          p2    p power.BSSR power.TRAD
# 1  0.2066667 0.006666667 0.14  0.9250322  0.9993007
# 2  0.2166667 0.016666667 0.15  0.8832410  0.9961241
# 3  0.2266667 0.026666667 0.16  0.8500245  0.9900890
# 4  0.2366667 0.036666667 0.17  0.8259856  0.9815361
# 5  0.2466667 0.046666667 0.18  0.8100250  0.9712275
# 6  0.2566667 0.056666667 0.19  0.8002297  0.9599332
# 7  0.2666667 0.066666667 0.20  0.7947331  0.9481556
# 8  0.2766667 0.076666667 0.21  0.7920642  0.9361171
# 9  0.2866667 0.086666667 0.22  0.7911416  0.9238884
# 10 0.2966667 0.096666667 0.23  0.7911888  0.9115032
# ...
########################################################################################################################
source('BinarySampleSize.R')
BinaryPowerBSSR = function(asmd.p1, asmd.p2, p, Delta.A, Delta.T, N1, N2, omega, r, alpha, tar.power, Test, restricted, weighted) {
  # Initial sample sizes used for BSSR
  N11 = ceiling(omega * N1)
  N12 = ceiling(omega * N2)
  # Estimate pooled proportion under blided fashion
  hat.p = outer(0:N11, 0:N12, '+') / (N11 + N12)
  # Unique values of \hat{p}
  U.hat.p = unique(c(hat.p))
  # Asign IDs which hat.p is matched to U.hat.p
  match.ID = matrix(match(hat.p, U.hat.p), nrow = N11 + 1, ncol = N12 + 1)
  # Derive hat{p}_{1} and hat{p}_{2} using hat{p} and Delta
  hat.p1 = pmin(1, U.hat.p + (1 / (1 + r)) * Delta.A)
  hat.p2 = pmax(0, U.hat.p - (r / (1 + r)) * Delta.A)
  # Sample size re-estimation
  BSSR.N2 = sapply(seq(U.hat.p), function(i) BinarySampleSize(hat.p1[i], hat.p2[i], r, alpha, tar.power, Test)[['N2']])
  if(restricted == TRUE) { N22 = pmax(N2, BSSR.N2) - N12 } else { N22 = pmax(N12, BSSR.N2) - N12 }
  N21 = r * N22
  # Final sample sizes
  hat.N1 = N11 + N21
  hat.N2 = N12 + N22
  hat.N = hat.N1 + hat.N2
  if(weighted == TRUE) {
    # Define weights
    prod.binom.prob = outer(dbinom(0:N11, N11, asmd.p1), dbinom(0:N12, N12, asmd.p2))
    w.h = as.numeric(tapply(prod.binom.prob, match.ID, FUN = sum))
    N.w = ceiling(sum(hat.N * w.h))
    # Final sample sizes
    hat.N = pmax(hat.N, N.w)
    hat.N2 = ceiling(hat.N / (1 + r))
    hat.N1 = ceiling(r * hat.N2) 
  }
  # True proportion for each treatment group
  p1 = p + (1 / (1 + r)) * Delta.T
  p2 = p - (r / (1 + r)) * Delta.T
  # Omit scenarios with negative or greater than 1 p_{j} values
  omit.ID = union(which(p1 %>>% 1), which(p2 %<<% 0))
  if(length(omit.ID) %!=% 0) {
    p1 = p1[-omit.ID]
    p2 = p2[-omit.ID]
    p = p[-omit.ID]
  }
  # Set probability mass functions of binomial distributions for the interim analysis
  dbinom1 = outer(X = 0:N11, Y = p1, function(X, Y) dbinom(X, N11, Y))
  dbinom2 = outer(X = 0:N12, Y = p2, function(X, Y) dbinom(X, N12, Y))
  # Powers given \hat{N}_{1} and \hat{N}_{2}
  power.stage2 = do.call(
    cbind,
    lapply(seq(U.hat.p), function(i) {
      # Set rejection region
      RR = BinaryRR(hat.N1[i], hat.N2[i], alpha, Test)
      # Return power
      sapply(seq(p), function(j) c(dbinom(0:hat.N1[i], hat.N1[i], p1[j]) %*% pbinom(rowSums(RR) - 1, hat.N2[i], p2[j])))
    })
  )
  # Final power for BSSR
  power.BSSR = sapply(seq(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * power.stage2[k, match(hat.p, U.hat.p)], na.rm = TRUE)
  })
  # Powers without BSSR (i.e., traditional design)
  power.TRAD = BinaryPower(p1, p2, N1, N2, alpha, Test)
  # Output
  out = data.frame(p1, p2, p, power.BSSR, power.TRAD)
  return(out)
}