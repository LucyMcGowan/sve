# Profile Likelihood Confidence Interval for SVE from Effect Measure (Single)

Computes profile likelihood CI for a single effect measure.

## Usage

``` r
sve_effect_profile_ci_single(theta_hat, var_log_theta, level)
```

## Arguments

- theta_hat:

  MLE of theta (scalar)

- var_log_theta:

  Variance of log(theta) (scalar)

- level:

  Confidence level

## Value

List with estimate, lower and upper confidence bounds
