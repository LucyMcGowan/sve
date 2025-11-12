#' Symmetric Vaccine Efficacy from Fitted Model
#'
#' Computes the symmetric vaccine efficacy (SVE) and associated confidence
#' intervals from a fitted regression model (Cox, Poisson, etc.) using
#' the exact likelihood from the model.
#'
#' @param model A fitted model object. Currently supports:
#'   * `coxph` objects from the survival package
#'   * `glm` objects with family poisson or binomial
#' @param effect_name Character. Name of the effect variable in the model.
#' @param data A data frame containing the original data used to fit the model.
#'   Required for profile likelihood method.
#' @param level Confidence level for the interval (default is 0.95).
#' @param method Method used to construct the confidence interval.
#'
#'   One of:
#'
#'   * `"profile"` (default): Profile likelihood confidence interval using the
#'     exact likelihood from the fitted model.
#'     Requires `data` argument.
#'
#'   * `"tanh-wald"`: Applies a hyperbolic arctangent
#'     transform before forming the Wald interval, then transforms back.
#'     Uses the variance extracted from the model.
#'
#'   * `"wald"`: Standard Wald interval on the untransformed scale.
#' @param smooth Logical. Should variance smoothing be applied near the null?
#'     Default is TRUE. Only used for Wald-based methods (`wald`, `tanh-wald`).
#'     Ignored for profile likelihood. Recommended to avoid instability when
#'     effect is near 1.
#' @param epsilon Numeric. The smoothing window. Only used for Wald-based
#'   methods (`wald`, `tanh-wald`) when `smooth = TRUE`. If `NULL` and
#'   `smooth = TRUE`, defaults to \eqn{z_{\alpha/2}
#'   \sqrt{\hat{p}_0(1-\hat{p}_0)/n_0 + \hat{p}_1(1-\hat{p}_1)/n_1}}
#'   where \eqn{\hat{p}_0 = x_0/n_0}, \eqn{\hat{p}_1 = x_1/n_1}, and
#'   \eqn{z_{\alpha/2}} is the critical value from the standard normal
#'   distribution corresponding to the confidence level.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{estimate}{The symmetric vaccine efficacy estimate.}
#'   \item{lower}{Lower bound of the confidence interval.}
#'   \item{upper}{Upper bound of the confidence interval.}
#'   \item{level}{Confidence interval level.}
#'   \item{method}{Indicates the method used for the confidence interval.}
#' }
#'
#' @details
#' This function extracts the coefficient and variance for the effect
#' variable from a fitted model and computes SVE. For the profile likelihood
#' method, it refits the model with the effect coefficient constrained
#' to values corresponding to different SVE values, using the exact likelihood
#' from the original model.
#'
#' The symmetric vaccine efficacy (SVE) from a relative effect measure is:
#' \deqn{
#' \text{SVE} = \frac{1 - \theta}{\max(1, \theta)}
#' }
#' where \eqn{\theta = \exp(\beta)} is the hazard ratio, relative risk, or
#' odds ratio from the model.
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' # Cox proportional hazards model
#' fit_cox <- coxph(Surv(time, status) ~ vaccination + age,
#'                  data = sim_trial_data)
#' sve_from_model(fit_cox, "vaccination", data = sim_trial_data)
#'
#' # Poisson regression
#' fit_poisson <- glm(status ~ vaccination + age,
#'                    family = poisson(), data = sim_trial_data)
#' sve_from_model(fit_poisson, "vaccination", data = sim_trial_data)
#'
#' # Wald methods don't require data
#' sve_from_model(fit_cox, "vaccination", method = "wald")
#' }
#'
#' @export
sve_from_model <- function(model,
                           effect_name,
                           data = NULL,
                           level = 0.95,
                           method = c("profile", "tanh-wald", "wald"),
                           smooth = TRUE,
                           epsilon = NULL) {
  method <- match.arg(method)
  check_confidence_level(level)

  # Validate effect_name exists in model
  coef_names <- names(stats::coef(model))
  if (!effect_name %in% coef_names) {
    cli::cli_abort(
      c("Effect variable {.var {effect_name}} not found in model.",
        "i" = "Available variables: {.var {coef_names}}")
    )
  }

  # Extract coefficient and variance
  beta_hat <- stats::coef(model)[effect_name]
  var_beta <- stats::vcov(model)[effect_name, effect_name]
  theta_hat <- exp(beta_hat)
  sve_hat <- sve_effect(theta_hat)

  # Profile method requires data
  if (method == "profile") {
    if (is.null(data)) {
      cli::cli_abort(
        c("Profile likelihood method requires the {.arg data} argument.",
          "i" = "Please provide the original data frame used to fit the model.",
          "i" = "Alternatively, use method = 'tanh-wald' or method = 'wald'.")
      )
    }

    # Validate that data contains the effect variable
    data_names <- names(data)
    if (!effect_name %in% data_names) {
      cli::cli_abort(
        c("Effect variable {.var {effect_name}} not found in provided data.",
          "i" = "Available variables: {.var {data_names}}")
      )
    }

    result <- sve_model_profile_ci(
      model = model,
      data = data,
      effect_name = effect_name,
      level = level
    )
  } else {
    # Wald methods don't need data
    result <- sve_from_effect(
      theta = theta_hat,
      var_log_theta = var_beta,
      method = method,
      level = level,
      smooth = smooth,
      epsilon = epsilon
    )
  }

  return(result)
}

#' Profile Likelihood CI for SVE from Model
#'
#' Computes profile likelihood confidence interval by refitting the model
#' with constrained effect coefficient.
#'
#' @param model Fitted model object
#' @param data Data frame used to fit the model
#' @param effect_name Name of effect variable
#' @param level Confidence level
#' @return Data frame with estimate and CI bounds
#' @keywords internal
sve_model_profile_ci <- function(model, data, effect_name, level) {
  beta_hat <- stats::coef(model)[effect_name]
  theta_hat <- exp(beta_hat)
  sve_hat <- sve_effect(theta_hat)

  ll_unconstrained <- as.numeric(stats::logLik(model))

  # Create constrained log-likelihood function
  constrained_ll_fn <- function(sve_target) {
    # Bound SVE away from extremes to avoid numerical issues
    sve_target <- pmax(pmin(sve_target, 0.9999), -0.9999)

    theta_s <- sve_to_theta(sve_target)
    if (is.na(theta_s) || theta_s <= 0 || !is.finite(theta_s)) {
      return(-Inf)
    }

    beta_s <- log(theta_s)
    if (!is.finite(beta_s)) {
      return(-Inf)
    }

    tryCatch({
      ll <- refit_model_constrained(model, data, effect_name, beta_s)
      return(ll)
    }, error = function(e) {
      return(-Inf)
    })
  }

  # Use generic profile likelihood CI function
  ci <- profile_likelihood_ci(
    param_hat = sve_hat,
    log_lik_unconstrained = ll_unconstrained,
    constrained_ll_fn = constrained_ll_fn,
    level = level,
    param_bounds = c(-0.9999, 0.9999),
    tol = 1e-8
  )

  data.frame(
    estimate = sve_hat,
    lower = ci$lower,
    upper = ci$upper,
    level = level,
    method = "Profile"
  )
}

#' Refit Model with Constrained Coefficient
#'
#' Refits a model with one coefficient held fixed at a specified value.
#'
#' @param model Original fitted model
#' @param data Data frame used to fit the original model
#' @param effect_name Name of variable to constrain
#' @param beta_fixed Value to fix the coefficient at
#' @return Log-likelihood of constrained model
#' @keywords internal
refit_model_constrained <- function(model, data, effect_name, beta_fixed) {
  UseMethod("refit_model_constrained")
}

#' @keywords internal
refit_model_constrained.coxph <- function(model, data, effect_name, beta_fixed) {
  formula_orig <- stats::formula(model)

  # Create working copy of data with offset
  data_work <- data
  data_work$.offset_var <- data_work[[effect_name]] * beta_fixed

  # Extract variable names from original formula
  terms_orig <- stats::terms(formula_orig)
  all_vars <- all.vars(stats::delete.response(terms_orig))
  other_vars <- setdiff(all_vars, effect_name)

  # Reconstruct response (handles Surv() objects correctly)
  response_call <- formula_orig[[2]]

  # Build new formula with offset
  if (length(other_vars) > 0) {
    rhs <- paste(c(other_vars, "offset(.offset_var)"), collapse = " + ")
  } else {
    rhs <- "offset(.offset_var)"
  }

  new_formula <- stats::reformulate(rhs, response = response_call)

  # Refit model
  fit_constrained <- survival::coxph(new_formula, data = data_work)

  as.numeric(stats::logLik(fit_constrained))
}

#' @keywords internal
refit_model_constrained.glm <- function(model, data, effect_name, beta_fixed) {
  formula_orig <- stats::formula(model)
  family_orig <- stats::family(model)

  # Calculate offset for constrained variable
  offset_new <- data[[effect_name]] * beta_fixed

  # Add to existing offset if present in original model
  if (!is.null(model$offset)) {
    offset_new <- offset_new + model$offset
  }

  # Extract predictor variables (excluding the constrained variable)
  terms_orig <- stats::terms(formula_orig)
  response <- formula_orig[[2]]
  all_vars <- all.vars(stats::delete.response(terms_orig))
  # Remove offset() calls and the effect variable
  predictor_vars <- setdiff(all_vars, c("offset", effect_name))

  # Build new formula
  if (length(predictor_vars) == 0) {
    new_formula <- stats::reformulate("1", response = response)
  } else {
    new_formula <- stats::reformulate(predictor_vars, response = response)
  }

  # Refit model with offset
  fit_constrained <- stats::glm(
    formula = new_formula,
    family = family_orig,
    data = data,
    offset = offset_new
  )

  as.numeric(stats::logLik(fit_constrained))
}

#' @keywords internal
refit_model_constrained.default <- function(model, data, effect_name, beta_fixed) {
  cli::cli_abort(
    c("Model class {.cls {class(model)}} is not supported.",
      "i" = "Supported classes: coxph, glm")
  )
}
