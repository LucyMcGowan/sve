# Profile Likelihood CI for SVE from Model

Computes profile likelihood confidence interval by refitting the model
with constrained effect coefficient.

## Usage

``` r
sve_model_profile_ci(model, data, effect_name, level)
```

## Arguments

- model:

  Fitted model object

- data:

  Data frame used to fit the model

- effect_name:

  Name of effect variable

- level:

  Confidence level

## Value

Data frame with estimate and CI bounds
