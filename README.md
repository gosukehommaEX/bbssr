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
library(tidyr)
library(ggplot2)

# Calculate bayesian decision probabilities
results = BayesDecisionProbBinary(
  prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0, gamma1 = 0.9, gamma2 = 0.3,
  pi1 = seq(0, 1, l = 101), pi2 = rep(0.2, 101), n1 = 12, n2 = 12, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL,
  m1 = 30, m2 = 30, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

# Display a figure showing bayesian decision making
figure = results %>% 
  as_tibble() %>% 
  mutate(theta = pi1 - pi2) %>% 
  pivot_longer(
    cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
  ) %>% 
  mutate(
    Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
  ) %>% 
  ggplot(aes(x = theta, y = Prob, group = Decision)) +
  geom_line(aes(color = Decision, linetype = Decision), linewidth = 2) +
  theme_bw() +
  scale_color_manual(
    values = c('Go' = '#658D1B', 'Gray' = '#939597', 'NoGo' = '#D91E49'),
    labels =  c('Go', 'Gray', 'NoGo')
  ) +
  scale_x_continuous(
    #expand=c(0, 0),
    limits = c(0 - 0.2, 1 - 0.2),
    breaks = seq(0 - 0.2, 1 - 0.2, l = 6)
  ) +
  scale_y_continuous(
    #expand=c(0, 0),
    limits = c(0, 1),
    breaks = seq(0, 1, l = 11)
  ) +
  labs(
    title = 'Probability of making Go/Gray/NoGo decision',
    x = expression(theta),
    y = 'Probability'
  ) +
  theme(
    text = element_text(size = 30),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.key.width = unit(4, 'cm'),
    legend.text = element_text(size = 30),
    legend.title = element_blank(), 
    legend.position = 'bottom'
  )
```

<img src="man/figures/README-figure.png" width="100%"/>

## References

Kang et al.(20XX)). Title
