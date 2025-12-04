# Symmetric Vaccine Efficacy from Relative Effect Measures

Computes the symmetric vaccine efficacy (SVE) and associated confidence
intervals from a relative effect measure (e.g., relative risk, hazard
ratio, odds ratio) obtained from a regression model.

## Usage

``` r
sve_from_effect(
  theta,
  var_log_theta,
  level = 0.95,
  method = c("profile", "tanh-wald", "wald")
)
```

## Arguments

- theta:

  Numeric. The relative effect measure (e.g., relative risk from Poisson
  regression, hazard ratio from Cox model, etc.). Can be a vector.

- var_log_theta:

  Numeric. The estimated variance of log(theta), typically obtained from
  model output (e.g., `vcov(model)` or the squared standard error). Must
  have the same length as `theta`.

- level:

  Confidence level for the interval (default is 0.95).

- method:

  Method used to construct the confidence interval.

  One of:

  - `"profile"` (default): Profile likelihood-based confidence interval
    using a normal approximation for log(theta). The interval inverts
    the likelihood ratio test after transforming to the SVE scale. To
    use the exact likelihood from a fitted regression model instead of
    the normal approximation, see
    [`sve_from_model()`](https://lucymcgowan.github.io/sve/reference/sve_from_model.md).

  - `"tanh-wald"`: Applies a hyperbolic arctangent transform before
    forming the Wald interval, then transforms back. Improves coverage
    when the estimate is near +/- 1.

  - `"wald"`: Standard Wald interval on the untransformed scale.

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

  Indicates the method used for the confidence interval.

## Details

The symmetric vaccine efficacy (SVE) from a relative effect measure is:
\$\$ \text{SVE} = \frac{1 - \theta}{\max(1, \theta)} \$\$ where
\\\theta\\ is the relative effect measure. When \\\theta \< 1\\, this
simplifies to \\\text{SVE} = 1 - \theta\\ (protective). When \\\theta \>
1\\, it becomes \\\text{SVE} = (1 - \theta)/\theta\\ (harmful).

## Examples

``` r
# Example with a hazard ratio from a Cox model
# Suppose: HR = 0.7, SE(log(HR)) = 0.15
hr <- 0.7
se_log_hr <- 0.15
var_log_hr <- se_log_hr^2
sve_from_effect(theta = hr, var_log_theta = var_log_hr)
#>   estimate      lower     upper level  method
#> 1      0.3 0.06075634 0.4783036  0.95 Profile

# With tanh-Wald method
sve_from_effect(theta = hr, var_log_theta = var_log_hr, method = "tanh-wald")
#>   estimate      lower     upper level    method
#> 1      0.3 0.08317729 0.4897028  0.95 tanh-Wald

# Example with multiple effect estimates (e.g., from subgroups)
hrs <- c(0.5, 0.8, 1.2)
var_log_hrs <- c(0.04, 0.025, 0.036)
sve_from_effect(theta = hrs, var_log_theta = var_log_hrs)
#>     estimate       lower     upper level  method
#> 1  0.5000000  0.26003647 0.6621455  0.95 Profile
#> 2  0.2000000 -0.08309734 0.4131823  0.95 Profile
#> 3 -0.1666667 -0.42546748 0.1726732  0.95 Profile

# Using output directly from a Cox model
if (FALSE) { # \dontrun{
library(survival)
fit <- coxph(Surv(time, status) ~ vaccination + age + baseline_risk,
             data = sim_trial_data)
hr <- exp(coef(fit)["vaccination"])
var_log_hr <- vcov(fit)["vaccination", "vaccination"]
sve_from_effect(theta = hr, var_log_theta = var_log_hr)
} # }
```
