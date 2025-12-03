#' Symmetric Vaccine Efficacy
#'
#' Computes the symmetric vaccine efficacy (SVE) and associated confidence
#' intervals based on observed proportions of events in vaccinated and
#' unvaccinated groups.
#'
#' By default, confidence intervals are constructed using the profile likelihood.
#' Two Wald-based intervals are available, one that uses the hyperbolic
#' arctangent (tanh-Wald) to keep the bounds within (-1, 1), and the
#' original-scale Wald. An exact method based on the Clopper-Pearson intervals
#' is also available.
#'
#' @param x0 Numeric. Number of events in the unvaccinated group.
#' @param x1 Numeric. Number of events in the vaccinated group.
#' @param n0 Integer. Sample size of the unvaccinated group.
#' @param n1 Integer. Sample size of the vaccinated group.
#' @param level Numeric. Confidence level for the interval (default is 0.95).
#' @param method Character. Method used to construct the confidence interval.
#'
#'   One of:
#'
#'   * `"profile"` (default): Profile likelihood confidence interval using likelihood
#'     ratio test inversion. More accurate than Wald-based methods, especially
#'     for small samples or extreme efficacy values.
#'
#'   * `"tanh-wald"`: Applies a hyperbolic arctangent
#'     transform before forming the Wald interval, then transforms back.
#'     Keeps the interval bounded in (-1,1).
#'
#'   * `"wald"`: Standard Wald interval on the untransformed scale.
#'
#'   * `"exact"`: Uses exact conditional inference by conditioning on total
#'     cases and inverting the conditional binomial distribution (Clopper-Pearson
#'     limits for the conditional parameter). Conservative due to discreteness,
#'     but reliable for small samples or boundary cases.
#' @param correction Logical. Whether to perform a bias correction, default: `TRUE`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{estimate}{The symmetric vaccine efficacy estimate.}
#'   \item{lower}{Lower bound of the confidence interval.}
#'   \item{upper}{Upper bound of the confidence interval.}
#'   \item{level}{Confidence interval level.}
#'   \item{method}{Indicates the method used.}
#' }
#'
#' @details
#' The symmetric vaccine efficacy (SVE) is defined as:
#' \deqn{
#' \text{SVE} = 2 \times \frac{p_0 - p_1}{p_0 + p_1 + |p_0 - p_1|}
#' }
#'
#' The profile likelihood method constructs confidence intervals by inverting
#' the likelihood ratio test. For a given SVE value, it finds the maximum
#' likelihood estimates of p0 and p1 subject to the constraint that their
#' SVE equals the specified value, then compares this constrained likelihood
#' to the unconstrained MLE.
#'
#' @examples
#' # Example: unvaccinated 10% infection rate, vaccinated 5%, equal group sizes
#'
#' # Profile likelihood method
#' est_sve(x0 = 10, x1 = 5, n0 = 100, n1 = 100, method = "profile")
#'
#' # Wald (tanh-Wald)
#' est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "tanh-wald")
#'
#' # Without transform (uses Wald interval)
#' est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "wald")
#'
#' # Exact method for small samples
#' est_sve(x0 = 10, x1 = 5, n0 = 100, n1 = 100, method = "exact")
#'
#'
#' @export
est_sve <- function(x0,
                    x1,
                    n0,
                    n1,
                    level = 0.95,
                    method = c("profile", "tanh-wald", "wald", "exact"),
                    correction = TRUE) {
  method <- match.arg(method)
  check_count_inputs(x0, x1, n0, n1)
  check_confidence_level(level)

  # Compute point estimates
  p0 <- x0 / n0
  p1 <- x1 / n1
  sve_val <- sve(p0, p1, n0, n1, correction)

  # Profile likelihood method
  if (method == "profile") {
    result <- sve_profile_ci(x0, x1, n0, n1, level, correction)
    return(result)
  }

  # Exact conditional method
  if (method == "exact") {
    alpha <- 1 - level
    m <- x0 + x1

    # Clopper-Pearson limits for conditional probability
    pi_L <- stats::qbeta(alpha / 2, x1, m - x1 + 1)
    pi_L[x1 == 0] <- 0

    pi_U <- stats::qbeta(1 - alpha / 2, x1 + 1, m - x1)
    pi_U[x1 == m] <- 1

    # Convert to relative risk scale
    theta_L <- (n0 / n1) * (pi_L / (1 - pi_L))
    theta_L[pi_L == 0] <- 0

    theta_U <- (n0 / n1) * (pi_U / (1 - pi_U))
    theta_U[pi_U == 1] <- Inf

    # Convert to SVE scale
    lower <- (1 - theta_U) / pmax(1, theta_U)
    lower[is.infinite(theta_U)] <- -1

    upper <- (1 - theta_L) / pmax(1, theta_L)
    upper[theta_L == 0] <- 1

    method_label <- "Exact"
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

    method_label <- "tanh-Wald"
  } else {
    var_sve <- sve_var(p0, p1, n0, n1)
    se_sve <- sqrt(var_sve)

    z_crit <- stats::qnorm(1 - (1 - level) / 2)
    lower <- sve_val - z_crit * se_sve
    upper <- sve_val + z_crit * se_sve

    method_label <- "Wald"
  }

  data.frame(
    estimate = sve_val,
    lower = lower,
    upper = upper,
    level = level,
    method = method_label
  )
}

#' Compute SVE from Proportions
#'
#' @param p0 Proportion in unvaccinated group
#' @param p1 Proportion in vaccinated group
#' @param n0 Sample size of unvaccinated group
#' @param n1 Sample size of vaccinated group
#' @param correction Whether to apply the bias correction
#' @return Symmetric vaccine efficacy
#' @keywords internal
sve <- function(p0, p1, n0, n1, correction) {
  check_proportions(p0, p1)
  sve <- (p0 - p1) / pmax(p0, p1)
  if (correction) {
    bias <- ifelse(p0 > p1,
                   -(p1 * (1 - p0)) / (n0 * p0^2),
                   (p0 * (1 - p1)) / (n1 * p1^2))
    sve <- sve - bias
  }
  return(sve)
}

#' Compute Variance of SVE via Delta Method
#'
#' @param p0 Proportion in unvaccinated group
#' @param p1 Proportion in vaccinated group
#' @param n0 Sample size of unvaccinated group
#' @param n1 Sample size of vaccinated group
#' @return Variance of SVE
#' @keywords internal
sve_var <- function(p0, p1, n0, n1) {
  check_proportions(p0, p1)

  sigma0 <- p0 * (1 - p0) / n0
  sigma1 <- p1 * (1 - p1) / n1

  result <- numeric(length(p0))
  abs_diff <- abs(p0 - p1)

  # Identify regions
  idx_p0_gt_p1 <- (p0 >= p1)
  idx_p1_gt_p0 <- (p1 > p0)

  # Case 1: p0 > p1
  if (any(idx_p0_gt_p1)) {
    result[idx_p0_gt_p1] <- (
      p1[idx_p0_gt_p1]^2 * sigma0[idx_p0_gt_p1] +
        p0[idx_p0_gt_p1]^2 * sigma1[idx_p0_gt_p1]
    ) / p0[idx_p0_gt_p1]^4
  }

  # Case 2: p1 > p0
  if (any(idx_p1_gt_p0)) {
    result[idx_p1_gt_p0] <- (
      p1[idx_p1_gt_p0]^2 * sigma0[idx_p1_gt_p0] +
        p0[idx_p1_gt_p0]^2 * sigma1[idx_p1_gt_p0]
    ) / p1[idx_p1_gt_p0]^4
  }
  result
}
