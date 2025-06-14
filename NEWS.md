# bbssr 1.0.2

## Minor Updates

* Fixed title case in DESCRIPTION file for CRAN submission
* Updated from "Re-estimation" to "Re-Estimation" as requested by CRAN

# bbssr 1.0.1

## Minor Updates

* **Function Removal**: Removed `ClopperPearsonCI()` function as it was not being used in the main BSSR functionality
* **Documentation Updates**: Updated all documentation to reflect the removal of confidence interval functionality
* **Package Optimization**: Streamlined package to focus on core BSSR methods

# bbssr 1.0.0

## Initial Release

This is the first release of `bbssr`, a comprehensive R package for blinded sample size re-estimation (BSSR) in two-arm clinical trials with binary endpoints.

### Main Features

* **Blinded Sample Size Re-estimation**: Implement adaptive trial designs with `BinaryPowerBSSR()`
* **Multiple Exact Statistical Tests**: Support for 5 different exact tests:
  - Pearson chi-squared test (`'Chisq'`)
  - Fisher exact test (`'Fisher'`)
  - Fisher mid-p test (`'Fisher-midP'`)
  - Z-pooled exact unconditional test (`'Z-pool'`)
  - Boschloo exact unconditional test (`'Boschloo'`)
* **Flexible Design Options**: Choose between restricted, unrestricted, and weighted BSSR approaches
* **Traditional Methods**: Calculate power (`BinaryPower()`) and sample sizes (`BinarySampleSize()`) for fixed designs
* **Exact Confidence Intervals**: Clopper-Pearson confidence intervals (`ClopperPearsonCI()`)
* **Rejection Regions**: Compute exact rejection regions (`BinaryRR()`)

### Design Approaches

* **Restricted Design**: Conservative approach ensuring final sample size ≥ initial sample size
* **Unrestricted Design**: Flexible approach allowing both sample size increases and decreases
* **Weighted Design**: Advanced approach using weighted averaging across interim scenarios

### Documentation

* Comprehensive documentation with examples for all functions
* Detailed vignettes explaining methodology and usage:
  - `vignette("bbssr-introduction")` - Getting started guide
  - `vignette("bbssr-statistical-methods")` - Statistical methodology
* Complete README with practical examples

### Statistical Validity

* All methods maintain exact Type I error control at specified α level
* Exact statistical tests rather than asymptotic approximations
* Suitable for small to moderate sample sizes common in clinical trials

### Dependencies

* Base R (≥ 3.5.0)
* `fpCompare` for robust floating-point comparisons
* `stats` for statistical functions

### Development

* Package follows R package development best practices
* Comprehensive documentation with roxygen2
* Ready for CRAN submission
