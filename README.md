
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sve

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
est_sve(p0 = 0.10, 
        p1 = 0.05, 
        n0 = 1000, 
        n1 = 1000)
#>   estimate     lower     upper level    method
#> 1      0.5 0.3191164 0.6457354  0.95 tanh-Wald
```

By default, confidence intervals are calculated on a transformed scale
and then back-transformed in order to ensure that they remain between -1
and 1. If you do not want to use this transformation, you can set the
option `transform = FALSE`:

``` r
est_sve(p0 = 0.10, 
        p1 = 0.05, 
        n0 = 1000, 
        n1 = 1000,
        transform = FALSE)
#>   estimate     lower     upper level method
#> 1      0.5 0.3360176 0.6639824  0.95   Wald
```

## Using adjusted effect measures

The `est_sve_adjusted()` function computes SVE from adjusted relative
effect measures (e.g., hazard ratios from Cox models, relative risks
from Poisson regression). Below is an example using a Cox proportional
hazards model.

``` r
library(survival)

fit <- coxph(Surv(time, status) ~ trt + celltype + age, data = veteran)
hr <- exp(coef(fit)["trt"])
var_log_hr <- vcov(fit)["trt", "trt"]

est_sve_adjusted(theta = hr, var_log_theta = var_log_hr)
#>       estimate      lower     upper level    method
#> trt -0.1639033 -0.5159536 0.2355141  0.95 tanh-Wald
```

## Methods overview

The symmetric vaccine efficacy (SVE) is defined as

$$\text{SVE}=\frac{2(p_0-p_1)}{p_0+p_1+|p_0-p_1|},$$

where:

- $p_0$ = event proportion in the unvaccinated group  
- $p_1$ = event proportion in the vaccinated group

This formulation ensures the estimator is bounded between -1 and 1.
Equivalently, it can be written in terms of a relative effect measure,
$\theta$ (such as a relative risk or hazard ratio):

$$\text{SVE} = \frac{2(1 - \theta)}{1 + \theta + |1 - \theta|}.$$
