# Profile Likelihood Confidence Interval for SVE from Effect Measure (Vectorized)

Computes profile likelihood CI for multiple effect measures.

## Usage

``` r
sve_effect_profile_ci(theta, var_log_theta, level)
```

## Arguments

- theta:

  MLE of theta (vector)

- var_log_theta:

  Variance of log(theta) (vector)

- level:

  Confidence level

## Value

Data frame with estimate and CI bounds
