# Symmetric Vaccine Efficacy from Fitted Model

Computes the symmetric vaccine efficacy (SVE) and associated confidence
intervals from a fitted regression model (Cox, Poisson, etc.) using the
exact likelihood from the model.

## Usage

``` r
sve_from_model(
  model,
  effect_name,
  data = NULL,
  level = 0.95,
  method = c("profile", "tanh-wald", "wald"),
  smooth = TRUE,
  epsilon = NULL
)
```

## Arguments

- model:

  A fitted model object. Currently supports:

  - `coxph` objects from the survival package

  - `glm` objects with family poisson or binomial

- effect_name:

  Character. Name of the effect variable in the model.

- data:

  A data frame containing the original data used to fit the model.
  Required for profile likelihood method.

- level:

  Confidence level for the interval (default is 0.95).

- method:

  Method used to construct the confidence interval.

  One of:

  - `"profile"` (default): Profile likelihood confidence interval using
    the exact likelihood from the fitted model. Requires `data`
    argument.

  - `"tanh-wald"`: Applies a hyperbolic arctangent transform before
    forming the Wald interval, then transforms back. Uses the variance
    extracted from the model.

  - `"wald"`: Standard Wald interval on the untransformed scale.

- smooth:

  Logical. Should variance smoothing be applied near the null? Default
  is TRUE. Only used for Wald-based methods (`wald`, `tanh-wald`).
  Ignored for profile likelihood. Recommended to avoid instability when
  effect is near 1.

- epsilon:

  Numeric. The smoothing window. Only used for Wald-based methods
  (`wald`, `tanh-wald`) when `smooth = TRUE`. If `NULL` and
  `smooth = TRUE`, defaults to \\z\_{\alpha/2}
  \sqrt{\hat{p}\_0(1-\hat{p}\_0)/n_0 + \hat{p}\_1(1-\hat{p}\_1)/n_1}\\
  where \\\hat{p}\_0 = x_0/n_0\\, \\\hat{p}\_1 = x_1/n_1\\, and
  \\z\_{\alpha/2}\\ is the critical value from the standard normal
  distribution corresponding to the confidence level.

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

This function extracts the coefficient and variance for the effect
variable from a fitted model and computes SVE. For the profile
likelihood method, it refits the model with the effect coefficient
constrained to values corresponding to different SVE values, using the
exact likelihood from the original model.

The symmetric vaccine efficacy (SVE) from a relative effect measure is:
\$\$ \text{SVE} = \frac{1 - \theta}{\max(1, \theta)} \$\$ where \\\theta
= \exp(\beta)\\ is the hazard ratio, relative risk, or odds ratio from
the model.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)

# Cox proportional hazards model
fit_cox <- coxph(Surv(time, status) ~ vaccination + age,
                 data = sim_trial_data)
sve_from_model(fit_cox, "vaccination", data = sim_trial_data)

# Poisson regression
fit_poisson <- glm(status ~ vaccination + age,
                   family = poisson(), data = sim_trial_data)
sve_from_model(fit_poisson, "vaccination", data = sim_trial_data)

# Wald methods don't require data
sve_from_model(fit_cox, "vaccination", method = "wald")
} # }
```
