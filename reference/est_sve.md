# Symmetric Vaccine Efficacy

Computes the symmetric vaccine efficacy (SVE) and associated confidence
intervals based on observed proportions of events in vaccinated and
unvaccinated groups.

## Usage

``` r
est_sve(
  x0,
  x1,
  n0,
  n1,
  level = 0.95,
  method = c("profile", "tanh-wald", "wald", "exact"),
  correction = FALSE
)
```

## Arguments

- x0:

  Numeric. Number of events in the unvaccinated group.

- x1:

  Numeric. Number of events in the vaccinated group.

- n0:

  Integer. Sample size of the unvaccinated group.

- n1:

  Integer. Sample size of the vaccinated group.

- level:

  Numeric. Confidence level for the interval (default is 0.95).

- method:

  Character. Method used to construct the confidence interval.

  One of:

  - `"profile"` (default): Profile likelihood confidence interval using
    likelihood ratio test inversion. More accurate than Wald-based
    methods, especially for small samples or extreme efficacy values.

  - `"tanh-wald"`: Applies a hyperbolic arctangent transform before
    forming the Wald interval, then transforms back. Keeps the interval
    bounded in (-1,1).

  - `"wald"`: Standard Wald interval on the untransformed scale.

  - `"exact"`: Uses exact conditional inference by conditioning on total
    cases and inverting the conditional binomial distribution
    (Clopper-Pearson limits for the conditional parameter). Conservative
    due to discreteness, but reliable for small samples or boundary
    cases.

- correction:

  Logical. Whether to perform a bias correction, default: `FALSE`.

## Value

A data frame with the following columns:

- estimate:

  The symmetric vaccine efficacy estimate.

- lower:

  Lower bound of the confidence interval.

- upper:

  Upper bound of the confidence interval.

- level:

  Confidence interval level.

- method:

  Indicates the method used.

## Details

By default, confidence intervals are constructed using the profile
likelihood. Two Wald-based intervals are available, one that uses the
hyperbolic arctangent (tanh-Wald) to keep the bounds within (-1, 1), and
the original-scale Wald. An exact method based on the Clopper-Pearson
intervals is also available.

The symmetric vaccine efficacy (SVE) is defined as: \$\$ \text{SVE} = 2
\times \frac{p_0 - p_1}{p_0 + p_1 + \|p_0 - p_1\|} \$\$

The profile likelihood method constructs confidence intervals by
inverting the likelihood ratio test. For a given SVE value, it finds the
maximum likelihood estimates of p0 and p1 subject to the constraint that
their SVE equals the specified value, then compares this constrained
likelihood to the unconstrained MLE.

## Examples

``` r
# Example: unvaccinated 10% infection rate, vaccinated 5%, equal group sizes

# Profile likelihood method
est_sve(x0 = 10, x1 = 5, n0 = 100, n1 = 100, method = "profile")
#>   estimate      lower     upper level  method
#> 1      0.5 -0.2614476 0.8394814  0.95 Profile

# Wald (tanh-Wald)
est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "tanh-wald")
#>   estimate     lower     upper level    method
#> 1      0.5 0.3191164 0.6457354  0.95 tanh-Wald

# Without transform (uses Wald interval)
est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "wald")
#>   estimate     lower     upper level method
#> 1      0.5 0.3360176 0.6639824  0.95   Wald

# Exact method for small samples
est_sve(x0 = 10, x1 = 5, n0 = 100, n1 = 100, method = "exact")
#>   estimate      lower     upper level method
#> 1      0.5 -0.3771404 0.8659031  0.95  Exact

```
