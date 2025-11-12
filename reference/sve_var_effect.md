# Compute Variance of SVE via Delta Method

Compute Variance of SVE via Delta Method

## Usage

``` r
sve_var_effect(theta, var_log_theta, smooth, epsilon, level)
```

## Arguments

- theta:

  Relative effect measure

- var_log_theta:

  Variance of log(theta)

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

Variance of SVE
