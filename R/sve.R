#' Symmetric Vaccine Efficacy
#'
#' Computes the **symmetric vaccine efficacy (SVE)** and associated confidence
#' intervals based on observed proportions of events in vaccinated and
#' unvaccinated groups.
#'
#' By default, confidence intervals are constructed using the hyperbolic
#' arctangent (tanh-Wald) for improved coverage, especially when efficacy values
#' approach the bounds of ±1.
#'
#' @param p0 Numeric. Proportion (risk or incidence) in the unvaccinated group.
#' @param p1 Numeric. Proportion (risk or incidence) in the vaccinated group.
#' @param n0 Integer. Sample size of the unvaccinated group.
#' @param n1 Integer. Sample size of the vaccinated group.
#' @param level Confidence level for the interval (default is 0.95).
#' @param transform Logical. If `TRUE` (default), applies a archtanh-transform
#'   before constructing confidence intervals and then transforms back to the
#'   original scale (tanh-Wald). If `FALSE`, uses a standard Wald interval.
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
#' The symmetric vaccine efficacy (SVE) is defined as:
#' \deqn{
#' \text{SVE} = 2 \times \frac{p_0 - p_1}{p_0 + p_1 + |p_0 - p_1|}
#' }
#' @examples
#' # Example: unvaccinated 10% infection rate, vaccinated 5%, equal group sizes
#' est_sve(p0 = 0.10, p1 = 0.05, n0 = 1000, n1 = 1000)
#'
#' # Without Fisher z-transform (uses Wald interval)
#' est_sve(p0 = 0.10, p1 = 0.05, n0 = 1000, n1 = 1000, transform = FALSE)
#'
#' @export
est_sve <- function(p0, p1, n0, n1, level = 0.95, transform = TRUE) {
  sve_val <- sve(p0, p1)

  if (transform) {
    z_val <- atanh(sve_val)

    var_sve <- sve_var(p0, p1, n0, n1)
    var_z <- var_sve / (1 - sve_val^2)^2
    se_z <- sqrt(var_z)

    z_crit <- stats::qnorm(1 - (1 - level) / 2)
    lower_z <- z_val - z_crit * se_z
    upper_z <- z_val + z_crit * se_z

    lower <- tanh(lower_z)
    upper <- tanh(upper_z)

  } else {
    var_sve <- sve_var(p0, p1, n0, n1)
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
    method = if (transform) "tanh-Wald" else "Wald"
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
