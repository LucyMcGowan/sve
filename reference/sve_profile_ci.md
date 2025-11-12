# Profile Likelihood Confidence Interval for SVE (Vectorized)

Computes profile likelihood CI for multiple sets of proportions.

## Usage

``` r
sve_profile_ci(x0, x1, n0, n1, level = 0.95)
```

## Arguments

- x0:

  Number of events in unvaccinated group (vector)

- x1:

  Number of events in vaccinated group (vector)

- n0:

  Sample size of unvaccinated group (vector)

- n1:

  Sample size of vaccinated group (vector)

- level:

  Confidence level

## Value

Data frame with estimate and CI bounds
