#' Profile Likelihood Confidence Interval for a Scalar Parameter
#'
#' Function to construct profile likelihood confidence intervals
#' by inverting the likelihood ratio test.
#'
#' @param param_hat MLE of the parameter of interest
#' @param log_lik_unconstrained Log-likelihood at the unconstrained MLE
#' @param constrained_ll_fn Function that takes a parameter value and returns
#'   the constrained (profile) log-likelihood at that value
#' @param level Confidence level (default 0.95)
#' @param param_bounds Two-element vector specifying search bounds for the
#'   parameter (default c(-0.9999, 0.9999))
#' @param tol Tolerance for root finding (default 1e-8)
#'
#' @return List with lower and upper confidence bounds
#' @keywords internal
profile_likelihood_ci <- function(param_hat,
                                  log_lik_unconstrained,
                                  constrained_ll_fn,
                                  level = 0.95,
                                  param_bounds = c(-0.9999, 0.9999),
                                  tol = 1e-8) {
  crit <- stats::qchisq(level, df = 1)

  lr_stat_fn <- function(param_val) {
    ll_constrained <- constrained_ll_fn(param_val)
    if (is.infinite(ll_constrained)) return(Inf)
    -2 * (ll_constrained - log_lik_unconstrained)
  }

  root_fn <- function(param_val) {
    lr_stat_fn(param_val) - crit
  }

  lower <- tryCatch({
    stats::uniroot(
      root_fn,
      interval = c(param_bounds[1], param_hat),
      tol = tol,
      extendInt = "downX"
    )$root
  }, error = function(e) {
    param_bounds[1]
  })

  upper <- tryCatch({
    stats::uniroot(
      root_fn,
      interval = c(param_hat, param_bounds[2]),
      tol = tol,
      extendInt = "upX"
    )$root
  }, error = function(e) {
    param_bounds[2]
  })

  list(lower = lower, upper = upper)
}

#' Binomial Log-Likelihood Helper
#'
#' Computes the binomial log-likelihood for a single proportion.
#'
#' @param x Number of successes
#' @param n Number of trials
#' @param p Probability of success
#' @return Log-likelihood value
#' @keywords internal
binomial_loglik <- function(x, n, p) {
  if (p <= 0 || p >= 1) return(-Inf)

  ll <- 0
  if (x > 0) ll <- ll + x * log(p)
  if (n - x > 0) ll <- ll + (n - x) * log(1 - p)

  ll
}

#' Two-Sample Binomial Log-Likelihood
#'
#' Computes the log-likelihood for two independent binomial samples.
#'
#' @param x0 Number of events in group 0
#' @param x1 Number of events in group 1
#' @param n0 Sample size of group 0
#' @param n1 Sample size of group 1
#' @param p0 Probability in group 0
#' @param p1 Probability in group 1
#'
#' @return Log-likelihood value
#' @keywords internal
two_sample_binomial_loglik <- function(x0, x1, n0, n1, p0, p1) {
  binomial_loglik(x0, n0, p0) + binomial_loglik(x1, n1, p1)
}

#' Convert SVE to Theta (Relative Risk)
#'
#' @param sve Symmetric vaccine efficacy value
#' @return Corresponding theta value
#' @keywords internal
sve_to_theta <- function(sve) {
  ifelse(sve >= 0, 1 - sve, 1 / (1 + sve))
}

#' Convert SVE to Proportions Given p0
#'
#' Given p0 and a target SVE, solve for p1.
#'
#' @param p0 Proportion in group 0
#' @param sve_target Target SVE value
#' @return Corresponding p1 value
#' @keywords internal
sve_to_p1 <- function(p0, sve_target) {
  if (sve_target >= 0) {
    p1 <- p0 * (1 - sve_target)
  } else {
    p1 <- p0 / (1 + sve_target)
  }

  p1
}

#' Determine Valid p0 Range for SVE Constraint
#'
#' @param sve_target Target SVE value
#' @param tol Tolerance for bounds (default 1e-8)
#' @return Two-element vector with lower and upper bounds for p0
#' @keywords internal
sve_p0_bounds <- function(sve_target, tol = 1e-8) {
  p0_lower <- tol
  p0_upper <- 1 - tol

  if (sve_target >= 0) {
    if (sve_target < 1) {
      p0_upper <- min(p0_upper, 1 / (1 - sve_target) - tol)
    }
  } else {
    p0_upper <- min(p0_upper, (1 + sve_target) - tol)
  }

  c(p0_lower, p0_upper)
}

#' Profile Likelihood Confidence Interval for SVE (Single Observation)
#'
#' Computes profile likelihood CI for a single set of proportions.
#'
#' @param x0 Number of events in unvaccinated group (scalar)
#' @param x1 Number of events in vaccinated group (scalar)
#' @param n0 Sample size of unvaccinated group (scalar)
#' @param n1 Sample size of vaccinated group (scalar)
#' @param level Confidence level
#' @param correction Logical. Whether to perform a bias correction, default: `TRUE`.
#' @return List with estimate and CI bounds
#' @keywords internal
sve_profile_ci_single <- function(x0, x1, n0, n1, level = 0.95,
                                  correction = TRUE) {
  # Handle edge cases where both groups have 0 or all events
  if (x0 == 0 && x1 == 0) {
    return(list(estimate = 0, lower = -1, upper = 1))
  }
  if (x0 == n0 && x1 == n1) {
    return(list(estimate = 0, lower = -1, upper = 1))
  }

  # Compute MLEs
  p0_hat <- x0 / n0
  p1_hat <- x1 / n1
  sve_hat <- sve(p0_hat, p1_hat, n0, n1, correction)

  # Unconstrained MLE log-likelihood
  ll_unconstrained <- two_sample_binomial_loglik(x0, x1, n0, n1, p0_hat, p1_hat)

  # Constrained log-likelihood as function of SVE
  constrained_ll_fn <- function(sve_target) {
    # Find optimal p0 given SVE constraint
    obj_fn <- function(p0) {
      p1 <- sve_to_p1(p0, sve_target)
      if (p1 <= 0 || p1 >= 1) return(-Inf)
      two_sample_binomial_loglik(x0, x1, n0, n1, p0, p1)
    }

    # Get valid p0 bounds for this SVE value
    p0_bounds <- sve_p0_bounds(sve_target, tol = 1e-8)
    if (p0_bounds[2] <= p0_bounds[1]) return(-Inf)

    opt <- tryCatch({
      stats::optimize(obj_fn, p0_bounds, maximum = TRUE, tol = 1e-10)
    }, error = function(e) {
      list(objective = -Inf)
    })

    opt$objective
  }

  ci <- profile_likelihood_ci(
    param_hat = sve_hat,
    log_lik_unconstrained = ll_unconstrained,
    constrained_ll_fn = constrained_ll_fn,
    level = level,
    param_bounds = c(-0.9999, 0.9999),
    tol = 1e-8
  )

  list(
    estimate = sve_hat,
    lower = ci$lower,
    upper = ci$upper
  )
}

#' Profile Likelihood Confidence Interval for SVE (Vectorized)
#'
#' Computes profile likelihood CI for multiple sets of proportions.
#'
#' @param x0 Number of events in unvaccinated group (vector)
#' @param x1 Number of events in vaccinated group (vector)
#' @param n0 Sample size of unvaccinated group (vector)
#' @param n1 Sample size of vaccinated group (vector)
#' @param level Confidence level
#' @param correction Logical. Whether to perform a bias correction, default: `TRUE`.
#' @return Data frame with estimate and CI bounds
#' @keywords internal
sve_profile_ci <- function(x0, x1, n0, n1, level = 0.95, correction = TRUE) {
  n <- length(x0)
  estimates <- numeric(n)
  lowers <- numeric(n)
  uppers <- numeric(n)

  for (i in seq_len(n)) {
    result <- sve_profile_ci_single(x0[i], x1[i], n0[i], n1[i], level, correction)
    estimates[i] <- result$estimate
    lowers[i] <- result$lower
    uppers[i] <- result$upper
  }

  data.frame(
    estimate = estimates,
    lower = lowers,
    upper = uppers,
    level = level,
    method = "Profile"
  )
}

#' Profile Likelihood Confidence Interval for SVE from Effect Measure (Single)
#'
#' Computes profile likelihood CI for a single effect measure.
#'
#' @param theta_hat MLE of theta (scalar)
#' @param var_log_theta Variance of log(theta) (scalar)
#' @param level Confidence level
#' @return List with estimate, lower and upper confidence bounds
#' @keywords internal
sve_effect_profile_ci_single <- function(theta_hat, var_log_theta, level) {
  # This assumes theta_hat follows a log-normal distribution
  log_lik_theta <- function(theta) {
    if (theta <= 0) return(-Inf)
    -0.5 * (log(theta) - log(theta_hat))^2 / var_log_theta
  }

  ll_unconstrained <- log_lik_theta(theta_hat)

  constrained_ll_fn <- function(sve_target) {
    theta_s <- sve_to_theta(sve_target)
    log_lik_theta(theta_s)
  }

  sve_hat <- sve_effect(theta_hat)

  ci <- profile_likelihood_ci(
    param_hat = sve_hat,
    log_lik_unconstrained = ll_unconstrained,
    constrained_ll_fn = constrained_ll_fn,
    level = level,
    param_bounds = c(-0.9999, 0.9999),
    tol = .Machine$double.eps^0.5
  )

  list(
    estimate = sve_hat,
    lower = ci$lower,
    upper = ci$upper
  )
}

#' Profile Likelihood Confidence Interval for SVE from Effect Measure (Vectorized)
#'
#' Computes profile likelihood CI for multiple effect measures.
#'
#' @param theta MLE of theta (vector)
#' @param var_log_theta Variance of log(theta) (vector)
#' @param level Confidence level
#' @return Data frame with estimate and CI bounds
#' @keywords internal
sve_effect_profile_ci <- function(theta, var_log_theta, level) {
  n <- length(theta)
  estimates <- numeric(n)
  lowers <- numeric(n)
  uppers <- numeric(n)

  for (i in seq_len(n)) {
    result <- sve_effect_profile_ci_single(theta[i], var_log_theta[i], level)
    estimates[i] <- result$estimate
    lowers[i] <- result$lower
    uppers[i] <- result$upper
  }

  data.frame(
    estimate = estimates,
    lower = lowers,
    upper = uppers,
    level = level,
    method = "Profile"
  )
}
