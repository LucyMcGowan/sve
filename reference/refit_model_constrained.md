# Refit Model with Constrained Coefficient

Refits a model with one coefficient held fixed at a specified value.

## Usage

``` r
refit_model_constrained(model, data, effect_name, beta_fixed)
```

## Arguments

- model:

  Original fitted model

- data:

  Data frame used to fit the original model

- effect_name:

  Name of variable to constrain

- beta_fixed:

  Value to fix the coefficient at

## Value

Log-likelihood of constrained model
