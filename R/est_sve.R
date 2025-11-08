#' Symmetric Vaccine Efficacy
#'
#' Computes the **symmetric vaccine efficacy (SVE)** and associated confidence
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
#' @param level Confidence level for the interval (default is 0.95).
#' @param method Method used to construct the confidence interval.
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
#'
#' @param c Correction parameter for variance smoothing (default 1.96).
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
est_sve <- function(x0, x1, n0, n1, level = 0.95,
                    method = c("profile", "tanh-wald", "wald", "exact"),
                    c = 1.96) {

  method <- match.arg(method)

  check_count_inputs(x_0, x_1, n_0, n_1)
  check_confidence_level(level)

  p0 <- x0 / n0
  p1 <- x1 / n1
  sve_val <- sve(p0, p1)

  if (method == "profile") {
    result <- sve_profile_likelihood_ci(x0, x1, n0, n1, level)
    return(result)
  }

  if (method == "exact") {
    alpha <- 1 - level
    m <- x0 + x1
    pi_L <- stats::qbeta(alpha / 2, x1, m - x1 + 1)
    pi_L[x1 == 0] <- 0
    pi_U <- stats::qbeta(1 - alpha / 2, x1 + 1, m - x1)
    pi_U[x1 == m] <- 1

    theta_L <- (n0 / n1) * (pi_L / (1 - pi_L))
    theta_L[pi_L == 0] <- 0
    theta_U <- (n0 / n1) * (pi_U / (1 - pi_U))
    theta_U[pi_U == 1] <- Inf

    lower <- (1 - theta_U) / pmax(1, theta_U)
    lower[is.infinite(theta_U)] <- -1
    upper <- (1 - theta_L) / pmax(1, theta_L)
    upper[theta_L == 0] <- 1

    method_name <- "Exact"
  } else if (method == "tanh-wald") {
    z_val <- atanh(sve_val)
    var_sve <- sve_var(p0, p1, n0, n1, c)
    var_z <- var_sve / (1 - sve_val^2)^2
    se_z <- sqrt(var_z)
    z_crit <- stats::qnorm(1 - (1 - level) / 2)

    lower_z <- z_val - z_crit * se_z
    upper_z <- z_val + z_crit * se_z
    lower <- tanh(lower_z)
    upper <- tanh(upper_z)

    method_name <- "tanh Wald"
  } else {
    var_sve <- sve_var(p0, p1, n0, n1, c)
    se_sve <- sqrt(var_sve)
    z_crit <- stats::qnorm(1 - (1 - level) / 2)

    lower <- sve_val - z_crit * se_sve
    upper <- sve_val + z_crit * se_sve

    method_name <- "Wald"
  }

  data.frame(
    estimate = sve_val,
    lower = lower,
    upper = upper,
    level = level,
    method = method_name
  )
}

#' Profile Likelihood Confidence Interval for SVE (Vectorized)
#'
#' @param x0 Number of events in unvaccinated group (vector)
#' @param x1 Number of events in vaccinated group (vector)
#' @param n0 Sample size of unvaccinated group (vector)
#' @param n1 Sample size of vaccinated group (vector)
#' @param level Confidence level
#' @return Data frame with estimate and CI bounds
#' @keywords internal
sve_profile_likelihood_ci <- function(x0, x1, n0, n1, level = 0.95) {

  n <- length(x0)
  estimates <- numeric(n)
  lowers <- numeric(n)
  uppers <- numeric(n)

  # Process each observation
  for (i in seq_len(n)) {
    result <- profile_likelihood_ci_single(x0[i], x1[i], n0[i], n1[i], level)
    estimates[i] <- result$estimate
    lowers[i] <- result$lower
    uppers[i] <- result$upper
  }

  data.frame(
    estimate = estimates,
    lower = lowers,
    upper = uppers,
    level = level,
    method = "Profile Likelihood"
  )
}

#' Profile Likelihood Confidence Interval for SVE (Single Observation)
#'
#' @param x0 Number of events in unvaccinated group (scalar)
#' @param x1 Number of events in vaccinated group (scalar)
#' @param n0 Sample size of unvaccinated group (scalar)
#' @param n1 Sample size of vaccinated group (scalar)
#' @param level Confidence level
#' @return List with estimate and CI bounds
#' @keywords internal
profile_likelihood_ci_single <- function(x0, x1, n0, n1, level = 0.95) {

  # Handle edge cases
  if (x0 == 0 && x1 == 0) {
    return(list(estimate = 0, lower = -1, upper = 1))
  }
  if (x0 == n0 && x1 == n1) {
    return(list(estimate = 0, lower = -1, upper = 1))
  }

  p0_hat <- x0 / n0
  p1_hat <- x1 / n1
  sve_hat <- sve(p0_hat, p1_hat)

  # Log-likelihood function
  loglik <- function(p0, p1) {
    if (p0 <= 0 || p0 >= 1 || p1 <= 0 || p1 >= 1) return(-Inf)

    ll <- 0
    if (x0 > 0) ll <- ll + x0 * log(p0)
    if (n0 - x0 > 0) ll <- ll + (n0 - x0) * log(1 - p0)
    if (x1 > 0) ll <- ll + x1 * log(p1)
    if (n1 - x1 > 0) ll <- ll + (n1 - x1) * log(1 - p1)

    return(ll)
  }

  # Unconstrained MLE log-likelihood
  ll_unconstrained <- loglik(p0_hat, p1_hat)

  # Solve for p1 given p0 and SVE target
  solve_for_p1 <- function(p0, sve_target) {
    if (sve_target >= 0) {
      p1 <- p0 * (1 - sve_target)
    } else {
      p1 <- p0 / (1 + sve_target)
    }
    return(p1)
  }

  # Constrained log-likelihood
  constrained_ll <- function(sve_target) {
    obj_fn <- function(p0) {
      p1 <- solve_for_p1(p0, sve_target)
      if (p1 <= 0 || p1 >= 1) return(-Inf)
      loglik(p0, p1)
    }

    # Determine valid range for p0
    if (sve_target >= 0) {
      p0_lower <- 1e-8
      p0_upper <- 1 - 1e-8

      if (sve_target < 1) {
        p0_upper <- min(p0_upper, 1 / (1 - sve_target) - 1e-8)
      }
    } else {
      p0_lower <- 1e-8
      p0_upper <- 1 - 1e-8

      p0_upper <- min(p0_upper, (1 + sve_target) - 1e-8)
    }

    if (p0_upper <= p0_lower) return(-Inf)

    opt <- tryCatch({
      stats::optimize(obj_fn, c(p0_lower, p0_upper), maximum = TRUE, tol = 1e-10)
    }, error = function(e) {
      list(objective = -Inf)
    })

    return(opt$objective)
  }

  # Likelihood ratio statistic
  lr_stat_fn <- function(sve_target) {
    ll_constrained <- constrained_ll(sve_target)
    if (is.infinite(ll_constrained)) return(Inf)
    -2 * (ll_constrained - ll_unconstrained)
  }

  crit <- stats::qchisq(level, df = 1)

  root_fn <- function(sve_target) {
    lr_stat_fn(sve_target) - crit
  }

  # Find lower bound
  lower <- tryCatch({
    stats::uniroot(root_fn, c(-0.9999, sve_hat), tol = 1e-8, extendInt = "downX")$root
  }, error = function(e) {
    -1
  })

  # Find upper bound
  upper <- tryCatch({
    stats::uniroot(root_fn, c(sve_hat, 0.9999), tol = 1e-8, extendInt = "upX")$root
  }, error = function(e) {
    1
  })

  list(estimate = sve_hat, lower = lower, upper = upper)
}

#' SVE calculation
#' @keywords internal
sve <- function(p0, p1) {
  check_proportions(p0, p1)
  2 * (p0 - p1) / (p0 + p1 + abs(p0 - p1))
}

#' SVE variance calculation
#' @keywords internal
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

  # |p0 - p1| <= epsilon (smoothing region)
  p_avg <- (p0[idx_smooth] + p1[idx_smooth]) / 2
  result[idx_smooth] <- (sigma0[idx_smooth] + sigma1[idx_smooth]) /
    (p_avg + epsilon[idx_smooth] / 4)^2

  return(result)
}
