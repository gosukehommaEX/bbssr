---
editor_options: 
  markdown: 
    wrap: 72
---

# bbssr: Blinded Sample Size Re-estimation for Binary Endpoints

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bbssr)](http://cran.r-project.org/package=bbssr)
[![R-CMD-check](https://github.com/gosukehommaEX/bbssr/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/bbssr/actions)
[![Codecov test
coverage](https://codecov.io/gh/gosukehommaEX/bbssr/branch/main/graph/badge.svg)](https://codecov.io/gh/gosukehommaEX/bbssr?branch=main)

<!-- badges: end -->

## Overview

`bbssr` is a comprehensive R package designed for **blinded sample size
re-estimation (BSSR)** in two-arm clinical trials with binary endpoints.
Unlike traditional fixed-sample designs, BSSR allows for adaptive sample
size adjustments during the trial while maintaining the statistical
integrity and blinding of the study.

### Key Features

-   **Blinded Sample Size Re-estimation**: Implement adaptive trial
    designs that adjust sample sizes based on pooled data without
    unblinding treatment assignments
-   **Multiple Exact Statistical Tests**: Support for five different
    exact statistical tests optimized for binary endpoints
-   **Flexible Design Options**: Choose between restricted,
    unrestricted, and weighted BSSR approaches
-   **Comprehensive Power Analysis**: Calculate exact power for both
    traditional and BSSR designs
-   **CRAN-Ready**: Fully documented with examples, vignettes, and
    comprehensive test coverage

### Statistical Methods Supported

The package implements five exact statistical tests specifically
designed for binary endpoints in clinical trials:

1.  **Pearson chi-squared test** (`'Chisq'`) - One-sided exact test
2.  **Fisher exact test** (`'Fisher'`) - Classical exact conditional
    test
3.  **Fisher mid-p test** (`'Fisher-midP'`) - Less conservative
    alternative to Fisher exact
4.  **Z-pooled exact unconditional test** (`'Z-pool'`) - Unconditional
    exact test with pooled variance
5.  **Boschloo exact unconditional test** (`'Boschloo'`) - Most powerful
    unconditional exact test

## Why Use BSSR?

Traditional clinical trials with fixed sample sizes often suffer from: -
**Inefficient resource allocation** when initial assumptions are
incorrect - **Underpowered studies** due to overly optimistic effect
size estimates - **Ethical concerns** about continuing underpowered or
overpowered trials

BSSR addresses these issues by: - **Maintaining statistical validity**
through exact methods - **Preserving blinding** by using only pooled
response rates - **Optimizing sample sizes** based on observed data -
**Improving trial efficiency** while controlling Type I error

## Installation

Install the released version from CRAN:

``` r
install.packages("bbssr")
```

Or install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/bbssr")
```

## Quick Start

### Basic Power Calculation

``` r
library(bbssr)

# Calculate power for a traditional design
power_traditional <- BinaryPower(
  p1 = 0.5,     # Response rate in treatment group
  p2 = 0.2,     # Response rate in control group  
  N1 = 40,      # Sample size in treatment group
  N2 = 40,      # Sample size in control group
  alpha = 0.025, # One-sided significance level
  Test = 'Fisher' # Statistical test
)

print(power_traditional)
```

### Sample Size Calculation

``` r
# Calculate required sample size
sample_size <- BinarySampleSize(
  p1 = 0.5,           # Expected response rate in treatment group
  p2 = 0.2,           # Expected response rate in control group
  r = 1,              # Allocation ratio (1:1)
  alpha = 0.025,      # One-sided significance level
  tar.power = 0.8,    # Target power
  Test = 'Boschloo'   # Most powerful exact test
)

print(sample_size)
```

### Blinded Sample Size Re-estimation

``` r
library(dplyr)

# BSSR with different design rules
bssr_result <- BinaryPowerBSSR(
  asmd.p1 = 0.45,      # Assumed response rate in treatment group
  asmd.p2 = 0.09,      # Assumed response rate in control group
  p = seq(0.1, 0.9, by = 0.1), # Range of pooled response rates
  Delta.A = 0.36,      # Assumed treatment effect
  Delta.T = 0.36,      # True treatment effect
  N1 = 24,             # Initial sample size in treatment group
  N2 = 24,             # Initial sample size in control group
  omega = 0.5,         # Fraction for interim analysis
  r = 1,               # Allocation ratio
  alpha = 0.025,       # Significance level
  tar.power = 0.8,     # Target power
  Test = 'Z-pool',     # Statistical test
  restricted = FALSE,   # Unrestricted design
  weighted = FALSE     # Non-weighted approach
)

head(bssr_result)
```

## Advanced Example: Comparing BSSR Designs

``` r
library(bbssr)
library(dplyr)
library(ggplot2)

# Compare different BSSR approaches
power_comparison <- tibble(
  Rule = factor(
    c('Restricted', 'Unrestricted', 'Weighted'), 
    levels = c('Restricted', 'Unrestricted', 'Weighted')
  ),
  restricted = c(TRUE, FALSE, FALSE),
  weighted = c(FALSE, FALSE, TRUE)
) %>% 
  group_by_all() %>% 
  reframe(r = c(1, 2), N1 = c(24, 36), N2 = c(24, 18)) %>% 
  group_by_all() %>% 
  reframe(
    BinaryPowerBSSR(
      asmd.p1 = 0.45, asmd.p2 = 0.09, p = seq(0.1, 0.9, by = 0.01),
      Delta.A = 0.36, Delta.T = 0.36, N1, N2, omega = 0.5, r, 
      alpha = 0.025, tar.power = 0.8, Test = 'Z-pool', 
      restricted, weighted
    )
  ) %>% 
  mutate(
    Rule = factor(Rule, levels = c('Restricted', 'Unrestricted', 'Weighted')),
    Allocation = paste0('Allocation ratio = ', r, ':1')
  )

# Visualize the results
ggplot(power_comparison, aes(x = p, y = power.BSSR, color = Rule)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Allocation) +
  geom_hline(yintercept = 0.8, color = 'gray', linetype = 'dashed') +
  labs(
    x = "Pooled Response Rate (θ)",
    y = "Power",
    title = "Power Comparison: BSSR Design Rules",
    subtitle = "Horizontal line shows target power = 0.8"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

## Key Functions

| Function             | Purpose                                              |
|--------------------------------------|----------------------------------|
| `BinaryPower()`      | Calculate power for traditional fixed-sample designs |
| `BinarySampleSize()` | Calculate required sample size for given power       |
| `BinaryPowerBSSR()`  | Calculate power for BSSR designs                     |
| `BinaryRR()`         | Compute rejection regions for exact tests            |
| `ClopperPearsonCI()` | Calculate exact confidence intervals                 |

## Design Options for BSSR

### Restricted Design (`restricted = TRUE`)

-   Conservative approach
-   Final sample size ≥ initially planned sample size
-   Maintains original study timeline

### Unrestricted Design (`restricted = FALSE`)

-   Flexible approach
-   Allows both increases and decreases in sample size
-   Optimizes efficiency based on observed data

### Weighted Design (`weighted = TRUE`)

-   Advanced approach
-   Uses weighted averaging across interim scenarios
-   Provides robust performance across different pooled response rates

## Statistical Validity

All methods in `bbssr` maintain exact Type I error control at the
specified α level. The package implements exact statistical tests rather
than asymptotic approximations, ensuring validity even for small sample
sizes commonly encountered in clinical trials.

## Vignettes

For detailed examples and theoretical background, see: -
`vignette("bbssr-introduction")` - Getting started guide -
`vignette("bbssr-statistical-methods")` - Statistical methodology -
`vignette("bbssr-design-comparison")` - Comparing BSSR approaches

## Citation

If you use `bbssr` in your research, please cite:

```         
Homma, G. (2025). bbssr: Blinded Sample Size Re-estimation for Binary Endpoints. 
R package version 1.0.0. https://CRAN.R-project.org/package=bbssr
```

## Contributing

Contributions are welcome! Please see our [Contributing
Guidelines](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the
[LICENSE.md](LICENSE.md) file for details.

## Author

**Gosuke Homma**

------------------------------------------------------------------------

**Note**: This package is designed for use by statisticians and clinical
researchers familiar with adaptive trial designs. For regulatory
submissions, please consult with biostatisticians and regulatory affairs
specialists to ensure compliance with relevant guidelines (FDA, EMA,
etc.).
