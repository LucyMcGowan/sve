
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sve

<!-- badges: start -->

[![R-CMD-check](https://github.com/LucyMcGowan/sve/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucyMcGowan/sve/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/LucyMcGowan/sve/graph/badge.svg)](https://app.codecov.io/gh/LucyMcGowan/sve)
<!-- badges: end -->

The goal of sve is to facilitate the calculation of symmetric vaccine
efficacy (SVE) and corresponding variance and confidence intervals.

## Installation

You can install the development version of sve like so:

``` r
devtools::install_github("LucyMcGowan/sve")
```

## Example

``` r
library(sve)
```

The core function, `est_sve()`, computes the symmetric vaccine efficacy
and its confidence interval given the event proportions and sample sizes
for the vaccinated and unvaccinated groups.

**Example:**

- unvaccinated 10% infection rate  
- vaccinated 5% infection rate  
- $n_0=n_1=1,000$

``` r
est_sve(x0 = 10, 
        x1 = 5, 
        n0 = 1000, 
        n1 = 1000)
#>   estimate      lower     upper level  method
#> 1      0.5 -0.2865842 0.8437684  0.95 Profile
```

By default, confidence intervals are calculated using the profile
likelihood method. To calculate Wald-type intervals on a transformed
scale and then back-transformed in order to ensure that they remain
between -1 and 1 use `method = "tanh-wald"`. If you do not want to use
this transformation, you can set the option `method = "wald"`:

``` r
est_sve(x0 = 10, 
        x1 = 5, 
        n0 = 1000, 
        n1 = 1000,
        method = "wald")
#>   estimate      lower    upper level method
#> 1      0.5 -0.3050449 1.305045  0.95   Wald
```

## Using relative effect measures

The `sve_from_model()` function computes SVE from relative effect
measures (e.g., hazard ratios from Cox models, relative risks from
Poisson regression) extracted from model objects. Below is an example
using a Cox proportional hazards model and a simulated data set provided
in this package (`sim_trial_data`).

``` r
library(survival)

fit <- coxph(Surv(time, status) ~ vaccination + age + baseline_risk, 
             data = sim_trial_data)

sve_from_model(model = fit, data = sim_trial_data, effect_name = "vaccination")
#>              estimate     lower     upper level  method
#> vaccination 0.6776701 0.5907393 0.7472103  0.95 Profile
```

## Methods overview

The symmetric vaccine efficacy (SVE) is defined as

$$\text{SVE}=\frac{(p_0-p_1)}{\max(p_0, p_1)},$$

where:

- $p_0$ = event proportion in the unvaccinated group  
- $p_1$ = event proportion in the vaccinated group

This formulation ensures the estimator is bounded between -1 and 1.
Equivalently, it can be written in terms of a relative effect measure,
$\theta$ (such as a relative risk or hazard ratio):

$$\text{SVE} = \frac{(1 - \theta)}{\max(1,\theta)}.$$
