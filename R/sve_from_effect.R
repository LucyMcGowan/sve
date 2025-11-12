#' Symmetric Vaccine Efficacy from Relative Effect Measures
#'
#' Computes the symmetric vaccine efficacy (SVE) and associated confidence
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
#'   * `"profile"` (default): Profile likelihood-based confidence interval
#'     using a normal approximation for log(theta). The interval inverts
#'     the likelihood ratio test after transforming to the SVE scale. To use
#'     the exact likelihood from a fitted regression model instead of the
#'     normal approximation, see [sve_from_model()].
#'
#'   * `"tanh-wald"`: Applies a hyperbolic arctangent
#'     transform before forming the Wald interval, then transforms back.
#'     Improves coverage when the estimate is near +/- 1.
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
#' The symmetric vaccine efficacy (SVE) from a relative effect measure is:
#' \deqn{
#' \text{SVE} = \frac{1 - \theta}{\max(1, \theta)}
#' }
#' where \eqn{\theta} is the relative effect measure. When \eqn{\theta < 1},
#' this simplifies to \eqn{\text{SVE} = 1 - \theta} (protective). When
#' \eqn{\theta > 1}, it becomes \eqn{\text{SVE} = (1 - \theta)/\theta} (harmful).
#'
#' @examples
#' # Example with a hazard ratio from a Cox model
#' # Suppose: HR = 0.7, SE(log(HR)) = 0.15
#' hr <- 0.7
#' se_log_hr <- 0.15
#' var_log_hr <- se_log_hr^2
#' sve_from_effect(theta = hr, var_log_theta = var_log_hr)
#'
#' # With tanh-Wald method
#' sve_from_effect(theta = hr, var_log_theta = var_log_hr, method = "tanh-wald")
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
sve_from_effect <- function(theta,
                            var_log_theta,
                            level = 0.95,
                            method = c("profile", "tanh-wald", "wald"),
                            smooth = TRUE,
                            epsilon = NULL) {
  method <- match.arg(method)
  check_confidence_level(level)
  check_theta(theta, var_log_theta)

  sve_val <- sve_effect(theta)

  if (method == "profile") {
    result <- sve_effect_profile_ci(theta, var_log_theta, level)
    lower <- result$lower
    upper <- result$upper
  } else {
    # Wald-based methods
    var_sve <- sve_var_effect(theta, var_log_theta, smooth, epsilon, level)

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
  }

  method_label <- switch(
    method,
    "tanh-wald" = "tanh-Wald",
    "wald" = "Wald",
    "profile" = "Profile"
  )

  data.frame(
    estimate = sve_val,
    lower = lower,
    upper = upper,
    level = level,
    method = method_label
  )
}

#' Compute SVE from Effect Measure
#'
#' @param theta Relative effect measure
#' @return Symmetric vaccine efficacy
#' @keywords internal
sve_effect <- function(theta) {
  (1 - theta) / pmax(1, theta)
}

#' Compute Variance of SVE via Delta Method
#'
#' @param theta Relative effect measure
#' @param var_log_theta Variance of log(theta)
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
#' @return Variance of SVE
#' @keywords internal
sve_var_effect <- function(theta, var_log_theta, smooth, epsilon, level) {
  if (smooth) {
    if (is.null(epsilon)) {
      c <- stats::qnorm(1 - (1 - level) / 2)
      epsilon <- c * sqrt(var_log_theta)
    }
  } else {
    if (!is.null(epsilon)) {
      cli::cli_warn(
        c("{.arg epsilon} is ignored when {.code smooth = FALSE}.",
          "i" = "Set {.code smooth = TRUE} to use the {.arg epsilon} parameter.")
      )
      epsilon <- 0
    }
    epsilon <- 0
  }
  result <- numeric(length(theta))

  abs_diff <- abs(1 - theta)

  # Identify regions
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
  if (any(idx_smooth)) {
    deriv_smooth <- -2 * theta[idx_smooth] /
      (1 + theta[idx_smooth] + epsilon[idx_smooth] / 2)
    result[idx_smooth] <- deriv_smooth^2 * var_log_theta[idx_smooth]
  }

  result
}
