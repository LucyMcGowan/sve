# Profile Likelihood Confidence Interval for a Scalar Parameter

Function to construct profile likelihood confidence intervals by
inverting the likelihood ratio test.

## Usage

``` r
profile_likelihood_ci(
  param_hat,
  log_lik_unconstrained,
  constrained_ll_fn,
  level = 0.95,
  param_bounds = c(-0.9999, 0.9999),
  tol = 1e-08
)
```

## Arguments

- param_hat:

  MLE of the parameter of interest

- log_lik_unconstrained:

  Log-likelihood at the unconstrained MLE

- constrained_ll_fn:

  Function that takes a parameter value and returns the constrained
  (profile) log-likelihood at that value

- level:

  Confidence level (default 0.95)

- param_bounds:

  Two-element vector specifying search bounds for the parameter (default
  c(-0.9999, 0.9999))

- tol:

  Tolerance for root finding (default 1e-8)

## Value

List with lower and upper confidence bounds
