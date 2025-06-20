% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BinaryPower.R
\name{BinaryPower}
\alias{BinaryPower}
\title{Power Calculation for Two-Arm Trials with Binary Endpoints}
\usage{
BinaryPower(p1, p2, N1, N2, alpha, Test)
}
\arguments{
\item{p1}{True probability of responders for group 1 (can be a vector with different values)}

\item{p2}{True probability of responders for group 2 (can be a vector with different values)}

\item{N1}{Sample size for group 1}

\item{N2}{Sample size for group 2}

\item{alpha}{One-sided level of significance}

\item{Test}{Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'}
}
\value{
A numeric value or vector of power values. If vectors are provided for p1 and p2,
a vector of powers corresponding to each combination will be returned.
}
\description{
Calculates power for two-arm trials with binary endpoints using exact statistical tests.
The function supports five different one-sided tests and can handle vectors of probabilities.
}
\details{
The function supports the following five one-sided tests:
\itemize{
\item The one-sided Pearson chi-squared test (Chisq)
\item The Fisher exact test (Fisher)
\item The Fisher mid-p test (Fisher-midP)
\item The Z-pooled exact unconditional test (Z-pool)
\item The Boschloo exact unconditional test (Boschloo)
}

The power calculation is based on the exact distribution of the test statistic
under the specified alternative hypothesis.
}
\examples{
# Simple power calculation with fast Chi-squared test
power1 <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 5, N2 = 5,
                     alpha = 0.025, Test = 'Chisq')
print(power1)

\donttest{
# More computationally intensive examples
# Single power calculation with larger sample size
power2 <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 10, N2 = 40,
                     alpha = 0.025, Test = 'Boschloo')
print(power2)

# Multiple power calculations
p1_vec <- c(0.5, 0.6, 0.7, 0.8)
p2_vec <- c(0.2, 0.2, 0.2, 0.2)
powers <- BinaryPower(p1 = p1_vec, p2 = p2_vec, N1 = 10, N2 = 40,
                     alpha = 0.025, Test = 'Fisher')
print(powers)
}

}
\author{
Gosuke Homma (\email{my.name.is.gosuke@gmail.com})
}
