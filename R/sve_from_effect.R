#' Symmetric Vaccine Efficacy from Relative Effect Measures
#'
#' Computes the **symmetric vaccine efficacy (SVE)** and associated confidence
#' intervals from a relative effect measure (e.g., relative risk,
#' hazard ratio, odds ratio) obtained from a regression model.
#'
#' @param theta Numeric. The relative effect measure (e.g., relative risk from
#'   Poisson regression, hazard ratio from Cox model, etc.). Can be a vector.
#' @param var_log_theta Numeric. The estimated variance of log(theta), typically
#'   obtained from model output (e.g., `vcov(model)` or the squared standard
#'   error). Must have the same length as `theta`.
#' @param level Confidence level for the interval (default is 0.95).
#' @param method Method used to construct the confidence interval.
#'
#'   One of:
#'
#'   * `"tanh-wald"` (default): Applies a hyperbolic arctangent
#'     transform before forming the Wald interval, then transforms back.
#'     Improves coverage when the estimate is near +/- 1.
#'
#'   * `"wald"`: Standard Wald interval on the untransformed scale.
#' @param c Numeric. Constant for determining epsilon in the smoothing
#'   approximation (default is 1.96). The smoothing parameter is set to
#'   `c * sqrt(var_log_theta)`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{estimate}{The symmetric vaccine efficacy estimate.}
#'   \item{lower}{Lower bound of the confidence interval.}
#'   \item{upper}{Upper bound of the confidence interval.}
#'   \item{level}{Confidence interval level.}
#'   \item{method}{Indicates whether the interval is a tanh-Wald interval or
#'   standard Wald.}
#' }
#'
#' @details
#' The symmetric vaccine efficacy (SVE) from a relative effect measure is:
#' \deqn{
#' \text{SVE} = \frac{2(1 - \theta)}{1 + \theta + |1 - \theta|}
#' }
#' where \eqn{\theta} is the relative effect measure.
#'
#' The variance is computed via the delta method applied to \eqn{\phi = \log(\theta)}.
#' When \eqn{|1 - \theta| > \varepsilon}, the derivative depends on whether
#' \eqn{\theta < 1} (protective) or \eqn{\theta > 1} (harmful). When
#' \eqn{|1 - \theta| \leq \varepsilon}, a smooth approximation is used.
#'
#' @examples
#' # Example with a hazard ratio from a Cox model
#' # Suppose: HR = 0.7, SE(log(HR)) = 0.15
#' hr <- 0.7
#' se_log_hr <- 0.15
#' var_log_hr <- se_log_hr^2
#' sve_from_effect(theta = hr, var_log_theta = var_log_hr)
#'
#' # Example with multiple effect estimates (e.g., from subgroups)
#' hrs <- c(0.5, 0.8, 1.2)
#' var_log_hrs <- c(0.04, 0.025, 0.036)
#' sve_from_effect(theta = hrs, var_log_theta = var_log_hrs)
#'
#' # Using output directly from a Cox model
#' \dontrun{
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ vaccination + age + baseline_risk,
#'              data = sim_trial_data)
#' hr <- exp(coef(fit)["vaccination"])
#' var_log_hr <- vcov(fit)["vaccination", "vaccination"]
#' sve_from_effect(theta = hr, var_log_theta = var_log_hr)
#' }
#'
#' @export
sve_from_effect <- function(theta, var_log_theta, level = 0.95,
                             method = c("tanh-wald", "wald"), c = 1.96) {
  method <- match.arg(method)
  if (any(theta <= 0)) {
    cli::cli_abort("{.arg theta} must be positive.")
  }
  if (any(var_log_theta < 0)) {
    cli::cli_abort("{.arg var_log_theta} must be non-negative.")
  }
  if (length(theta) != length(var_log_theta)) {
    cli::cli_abort(
      c("{.arg theta} and {.arg var_log_theta} must have the same length.",
        "i" = "{.arg theta} has length {length(theta)}.",
        "i" = "{.arg var_log_theta} has length {length(var_log_theta)}.")
    )
  }

  sve_val <- sve_effect(theta)

  var_sve <- sve_var_effect(theta, var_log_theta, c = c)

  if (method == "tanh-wald") {
    z_val <- atanh(sve_val)
    var_z <- var_sve / (1 - sve_val^2)^2
    se_z <- sqrt(var_z)
    z_crit <- stats::qnorm(1 - (1 - level) / 2)
    lower_z <- z_val - z_crit * se_z
    upper_z <- z_val + z_crit * se_z
    lower <- tanh(lower_z)
    upper <- tanh(upper_z)
  } else {
    se_sve <- sqrt(var_sve)
    z_crit <- stats::qnorm(1 - (1 - level) / 2)
    lower <- sve_val - z_crit * se_sve
    upper <- sve_val + z_crit * se_sve
  }

  data.frame(
    estimate = sve_val,
    lower = lower,
    upper = upper,
    level = level,
    method = if (method == "tanh-wald") "tanh-Wald" else "Wald"
  )
}

sve_effect <- function(theta) {
  2 * (1 - theta) / (1 + theta + abs(1 - theta))
}

sve_var_effect <- function(theta, var_log_theta, c = 1.96) {
  epsilon <- c * sqrt(var_log_theta)

  result <- numeric(length(theta))
  abs_diff <- abs(1 - theta)

  idx_smooth <- abs_diff <= epsilon
  idx_protective <- (theta < 1) & !idx_smooth
  idx_harmful <- (theta > 1) & !idx_smooth

  # Case 1: theta < 1 (protective) and |1 - theta| > epsilon
  # d(SVE)/d(phi) = -theta
  result[idx_protective] <- theta[idx_protective]^2 * var_log_theta[idx_protective]

  # Case 2: theta > 1 (harmful) and |1 - theta| > epsilon
  # d(SVE)/d(phi) = -1/theta
  result[idx_harmful] <- (1 / theta[idx_harmful])^2 * var_log_theta[idx_harmful]

  # Case 3: |1 - theta| <= epsilon (near null)
  # d(SVE)/d(phi) = -2*theta / (1 + theta + epsilon/2)
  deriv_smooth <- -2 * theta[idx_smooth] / (1 + theta[idx_smooth] + epsilon[idx_smooth]/2)
  result[idx_smooth] <- deriv_smooth^2 * var_log_theta[idx_smooth]

  return(result)
}
