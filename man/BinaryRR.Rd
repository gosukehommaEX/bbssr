% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BinaryRR.R
\name{BinaryRR}
\alias{BinaryRR}
\title{Rejection Region for Two-Arm Trials with Binary Endpoints}
\usage{
BinaryRR(N1, N2, alpha, Test)
}
\arguments{
\item{N1}{Sample size for group 1}

\item{N2}{Sample size for group 2}

\item{alpha}{One-sided level of significance}

\item{Test}{Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'}
}
\value{
A logical matrix representing the rejection region (RR). Matrix dimensions
are (N1+1) x (N2+1), where TRUE indicates rejection of the null hypothesis.
}
\description{
Provides a rejection region (RR) for two-arm trials with binary endpoints using
various exact statistical tests. The function supports five different one-sided tests.
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
}
\examples{
# Simple example with small sample sizes (runs quickly)
N1 <- 5
N2 <- 5
alpha <- 0.025
Test <- 'Chisq'
RR <- BinaryRR(N1, N2, alpha, Test)
print(dim(RR))  # Should be (6, 6)

\donttest{
# More computationally intensive example
N1 <- 20
N2 <- 10
alpha <- 0.025
Test <- 'Boschloo'
RR <- BinaryRR(N1, N2, alpha, Test)
print(RR)
}

}
\author{
Gosuke Homma (\email{my.name.is.gosuke@gmail.com})
}
