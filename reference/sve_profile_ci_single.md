# Profile Likelihood Confidence Interval for SVE (Single Observation)

Computes profile likelihood CI for a single set of proportions.

## Usage

``` r
sve_profile_ci_single(x0, x1, n0, n1, level = 0.95)
```

## Arguments

- x0:

  Number of events in unvaccinated group (scalar)

- x1:

  Number of events in vaccinated group (scalar)

- n0:

  Sample size of unvaccinated group (scalar)

- n1:

  Sample size of vaccinated group (scalar)

- level:

  Confidence level

## Value

List with estimate and CI bounds
