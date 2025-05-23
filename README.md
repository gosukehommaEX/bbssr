# ExactBinaryFamily

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ExactBinaryFamily)](http://cran.r-project.org/package=ExactBinaryFamily)
<!-- badges: end -->

## Overview

`ExactBinaryFamily` is an R package that provides methods to Provides functions for calculating power, sample size, and rejection regions
    for two-arm trials with binary endpoints using exact statistical tests. Supports
    five different test methods: Pearson chi-squared test, Fisher exact test, 
    Fisher mid-p test, Z-pooled exact unconditional test, and Boschloo exact 
    unconditional test. Also includes functionality for blinded sample size 
    re-estimation (BSSR).

For technical details about the methodology, please refer to xxx et al.(20XX).

## Installation

You can install the development version of ExactBinaryFamily from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/ExactBinaryFamily")
```

## Usage

``` r
# Load packages
library(ExactBinaryFamily)
library(dplyr)
library(ggplot2)

# Calculate power of the BSSR with binary endpoint
power.BSSR = tibble(
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
      asmd.p1 = 0.45, asmd.p2 = 0.09, p = seq(0, 1, by = 0.01),
      Delta.A = 0.36, Delta.T = 0.36, N1, N2, omega = 0.5, r, 
      alpha = 0.025, tar.power = 0.8, Test = 'Z-pool', restricted, weighted
    )
  ) %>% 
  mutate(
    Rule = factor(Rule, levels = c('Restricted', 'Unrestricted', 'Weighted')),
    r = paste0('r==', r)
  ) %>% 
  rename(Theta = p, Power = power.BSSR)

# Display a figure
power.BSSR %>% 
  ggplot(aes(x = Theta, y = Power)) +
  geom_line(aes(color = Rule), linewidth = 1.2) +
  theme_bw() +
  facet_grid(. ~ r, labeller = label_parsed) + 
  geom_hline(yintercept = 0.8, color = 'gray', linetype = 'longdash', linewidth = 1.2) +
  scale_x_continuous(
    limits = c(0.1, 0.9),
    breaks = seq(0.1, 0.9, l = 5)
  ) +
  scale_y_continuous(
    limits = c(0.7, 1.0),
    breaks = seq(0.7, 1.0, l = 7)
  ) +
  labs(x = expression(theta)) +
  theme(
    text = element_text(size = 40),
    panel.spacing = unit(-0.1, 'lines'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 40),
    legend.title = element_blank(),
    legend.title.position = 'top',
    legend.position = 'bottom'
  )
```

<img src="man/figures/fig.README.png" width="100%"/>

## References

xxx et al.(20XX). Title
