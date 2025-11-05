#' Symmetric Vaccine Efficacy
#'
#' Computes the **symmetric vaccine efficacy (SVE)** and associated confidence
#' intervals based on observed proportions of events in vaccinated and
#' unvaccinated groups.
#'
#' By default, confidence intervals are constructed using the hyperbolic
#' arctangent (tanh-Wald) for improved coverage, especially when efficacy values
#' approach the bounds of +/- 1. An exact method based on the Clopper-Pearson
#' intervals is also available for small samples or boundary cases.
#'
#' @param x0 Numeric. Number of events in the unvaccinated group.
#' @param x1 Numeric. Number of events in the vaccinated group.
#' @param n0 Integer. Sample size of the unvaccinated group.
#' @param n1 Integer. Sample size of the vaccinated group.
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
#'
#'   * `"exact"`: Uses Clopper–Pearson binomial intervals for each group and
#'     propagates them through the SVE function. Conservative, but reliable
#'     for small samples or boundary cases.
#' @return A data frame with the following columns:
#' \describe{
#'   \item{estimate}{The symmetric vaccine efficacy estimate.}
#'   \item{lower}{Lower bound of the confidence interval.}
#'   \item{upper}{Upper bound of the confidence interval.}
#'   \item{level}{Confidence interval level.}
#'   \item{method}{Indicates the method used: "tanh-Wald", "Wald", or "Exact".}
#' }
#'
#' @details
#' The symmetric vaccine efficacy (SVE) is defined as:
#' \deqn{
#' \text{SVE} = 2 \times \frac{p_0 - p_1}{p_0 + p_1 + |p_0 - p_1|}
#' }
#' @examples
#' # Example: unvaccinated 10% infection rate, vaccinated 5%, equal group sizes
#' est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "tanh-wald")
#'
#' # Without Fisher z-transform (uses Wald interval)
#' est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "wald")
#'
#' # Exact method for small samples
#' est_sve(x0 = 10, x1 = 5, n0 = 100, n1 = 100, method = "exact")
#' @export
est_sve <- function(x0, x1, n0, n1, level = 0.95, method = c("tanh-wald", "wald", "exact")) {

  method <- match.arg(method)

  p0 <- x0 / n0
  p1 <- x1 / n1
  sve_val <- sve(p0, p1)

  if (method == "exact") {

    alpha <- 1 - level
    m <- x0 + x1

    if (x1 == 0) {
      pi_L <- 0
    } else {
      pi_L <- stats::qbeta(alpha / 2, x1, m - x1 + 1)
    }
    if (x1 == m) {
      pi_U <- 1
    } else {
      pi_U <- stats::qbeta(1 - alpha / 2, x1 + 1, m - x1)
    }

    if (pi_L == 0) {
      theta_L <- 0
    } else {
      theta_L <- (n0 / n1) * (pi_L / (1 - pi_L))
    }
    if (pi_U == 1) {
      theta_U <- Inf
    } else {
      theta_U <- (n0 / n1) * (pi_U / (1 - pi_U))
    }

    if (is.infinite(theta_U)) {
      lower <- -1
    } else {
      lower <- (1 - theta_U) / max(1, theta_U)
    }
    if (theta_L == 0) {
      upper <- 1
    } else {
      upper <- (1 - theta_L) / max(1, theta_L)
    }

    method <- "Exact"

  } else if (method == "tanh-wald") {
    z_val <- atanh(sve_val)

    var_sve <- sve_var(p0, p1, n0, n1)
    var_z <- var_sve / (1 - sve_val^2)^2
    se_z <- sqrt(var_z)

    z_crit <- stats::qnorm(1 - (1 - level) / 2)
    lower_z <- z_val - z_crit * se_z
    upper_z <- z_val + z_crit * se_z

    lower <- tanh(lower_z)
    upper <- tanh(upper_z)

    method <- "tanh Wald"

  } else {
    var_sve <- sve_var(p0, p1, n0, n1)
    se_sve <- sqrt(var_sve)
    z_crit <- stats::qnorm(1 - (1 - level) / 2)

    lower <- sve_val - z_crit * se_sve
    upper <- sve_val + z_crit * se_sve
    method <- "Wald"
  }

  data.frame(
    estimate = sve_val,
    lower = lower,
    upper = upper,
    level = level,
    method = method
  )
}

sve <- function(p0, p1) {
  check_proportions(p0, p1)
  2 * (p0 - p1) / (p0 + p1 + abs(p0 - p1))
}

sve_var <- function(p0, p1, n0, n1, c = 1.96) {
  check_proportions(p0, p1)
  sigma0 <- p0 * (1 - p0) / n0
  sigma1 <- p1 * (1 - p1) / n1

  epsilon <- c * sqrt(sigma0 + sigma1)
  result <- numeric(length(p0))

  abs_diff <- abs(p0 - p1)
  idx_smooth <- abs_diff <= epsilon
  idx_p0_gt_p1 <- (p0 > p1) & !idx_smooth
  idx_p1_gt_p0 <- (p1 > p0) & !idx_smooth

  # p0 > p1 and |p0 - p1| > epsilon
  result[idx_p0_gt_p1] <- (p1[idx_p0_gt_p1]^2 * sigma0[idx_p0_gt_p1] +
                             p0[idx_p0_gt_p1]^2 * sigma1[idx_p0_gt_p1]) /
    p0[idx_p0_gt_p1]^4

  # p1 > p0 and |p0 - p1| > epsilon
  result[idx_p1_gt_p0] <- (p1[idx_p1_gt_p0]^2 * sigma0[idx_p1_gt_p0] +
                             p0[idx_p1_gt_p0]^2 * sigma1[idx_p1_gt_p0]) /
    p1[idx_p1_gt_p0]^4

  # (|p0 - p1| <= epsilon)

  p_avg <- (p0[idx_smooth] + p1[idx_smooth]) / 2
  result[idx_smooth] <- (sigma0[idx_smooth] + sigma1[idx_smooth]) /
    (p_avg + epsilon[idx_smooth]/4)^2

  return(result)
}
