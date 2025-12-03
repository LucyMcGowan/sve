test_that("profile_likelihood_ci validates and handles basic cases", {
  # Simple quadratic log-likelihood
  constrained_ll <- function(x) -0.5 * (x - 2)^2

  result <- profile_likelihood_ci(
    param_hat = 2,
    log_lik_unconstrained = constrained_ll(2),
    constrained_ll_fn = constrained_ll,
    level = 0.95,
    param_bounds = c(-5, 5)
  )

  expect_true(result$lower < 2)
  expect_true(result$upper > 2)
  expect_true(is.finite(result$lower))
  expect_true(is.finite(result$upper))
})

test_that("sve_to_theta correctly transforms SVE", {
  # Protective: SVE >= 0
  expect_equal(sve_to_theta(0), 1)
  expect_equal(sve_to_theta(0.5), 0.5)
  expect_equal(sve_to_theta(0.9), 0.1)

  # Harmful: SVE < 0
  expect_equal(sve_to_theta(-0.5), 2)
  expect_equal(sve_to_theta(-0.99), 100)
})

test_that("sve_to_p1 correctly computes p1 from p0 and SVE", {
  # Protective effect
  expect_equal(sve_to_p1(0.1, 0.5), 0.05)
  expect_equal(sve_to_p1(0.2, 0.5), 0.1)

  # Harmful effect
  expect_equal(sve_to_p1(0.1, -0.5), 0.2)

  # Null effect
  expect_equal(sve_to_p1(0.1, 0), 0.1)
})

test_that("sve_profile_ci handles moderate efficacy", {
  result <- sve_profile_ci(
    x0 = 100, x1 = 50, n0 = 1000, n1 = 1000,
    level = 0.95,
    correction = FALSE
  )

  expect_equal(result$estimate, 0.5)
  expect_equal(result$method, "Profile")
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
  expect_gt(result$lower, -1)
  expect_lt(result$upper, 1)
})

test_that("sve_profile_ci handles high efficacy", {
  result <- sve_profile_ci(
    x0 = 100, x1 = 5, n0 = 1000, n1 = 1000,
    level = 0.95,
    correction = FALSE
  )

  expect_gt(result$estimate, 0.9)
  expect_true(is.finite(result$lower))
  expect_true(is.finite(result$upper))
  expect_gt(result$lower, 0.8)
})

test_that("sve_profile_ci handles harmful effects", {
  result <- sve_profile_ci(
    x0 = 50, x1 = 100, n0 = 1000, n1 = 1000,
    level = 0.95,
    correction = FALSE
  )

  expect_equal(result$estimate, -0.5)
  expect_lt(result$estimate, 0)
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("sve_profile_ci handles edge case: both groups have no events", {
  result <- sve_profile_ci(
    x0 = 0, x1 = 0, n0 = 100, n1 = 100,
    level = 0.95
  )

  expect_equal(result$estimate, 0)
  expect_equal(result$lower, -1)
  expect_equal(result$upper, 1)
})

test_that("sve_profile_ci handles edge case: both groups have all events", {
  result <- sve_profile_ci(
    x0 = 100, x1 = 100, n0 = 100, n1 = 100,
    level = 0.95
  )

  expect_equal(result$estimate, 0)
  expect_equal(result$lower, -1)
  expect_equal(result$upper, 1)
})

test_that("sve_profile_ci is vectorized", {
  result <- sve_profile_ci(
    x0 = c(100, 50, 30),
    x1 = c(50, 100, 60),
    n0 = c(1000, 1000, 1000),
    n1 = c(1000, 1000, 1000),
    level = 0.95
  )

  expect_equal(nrow(result), 3)
  expect_gt(result$estimate[1], 0)  # Protective
  expect_lt(result$estimate[2], 0)  # Harmful
  expect_lt(result$estimate[3], 0)  # Harmful
})

test_that("sve_effect_profile_ci handles moderate effects", {
  result <- sve_effect_profile_ci(
    theta = 0.7,
    var_log_theta = 0.0225,
    level = 0.95
  )

  expect_equal(result$estimate, 0.3, tolerance = 1e-10)
  expect_equal(result$method, "Profile")
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("sve_effect_profile_ci handles effects near null", {
  result <- sve_effect_profile_ci(
    theta = 0.95,
    var_log_theta = 0.01,
    level = 0.95
  )

  expect_equal(result$estimate, 0.05, tolerance = 1e-10)
  expect_true(is.finite(result$lower))
  expect_true(is.finite(result$upper))
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("sve_effect_profile_ci handles harmful effects", {
  result <- sve_effect_profile_ci(
    theta = 1.5,
    var_log_theta = 0.04,
    level = 0.95
  )

  expect_equal(result$estimate, -1/3, tolerance = 1e-10)
  expect_lt(result$estimate, 0)
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("sve_effect_profile_ci is vectorized", {
  result <- sve_effect_profile_ci(
    theta = c(0.5, 0.8, 1.2),
    var_log_theta = c(0.04, 0.025, 0.036),
    level = 0.95
  )

  expect_equal(nrow(result), 3)
  expect_gt(result$estimate[1], 0)  # Protective
  expect_gt(result$estimate[2], 0)  # Protective
  expect_lt(result$estimate[3], 0)  # Harmful
})

# Tests for sve_from_model.R

test_that("sve_from_model validates effect_name exists", {
  skip_if_not_installed("survival")

  data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    status = c(1, 1, 0, 1, 0, 1, 1, 0, 1, 0),
    treatment = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0),
    age = c(50, 55, 60, 65, 70, 52, 58, 62, 68, 72)
  )

  fit <- survival::coxph(survival::Surv(time, status) ~ treatment + age, data = data)

  expect_error(
    sve_from_model(fit, "vaccine", data = data),
    "not found in model"
  )
})

test_that("sve_from_model requires data for profile method", {
  skip_if_not_installed("survival")

  data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    status = c(1, 1, 0, 1, 0, 1, 1, 0, 1, 0),
    treatment = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0),
    age = c(50, 55, 60, 65, 70, 52, 58, 62, 68, 72)
  )

  fit <- survival::coxph(survival::Surv(time, status) ~ treatment + age, data = data)

  expect_error(
    sve_from_model(fit, "treatment", method = "profile"),
    "requires.*data"
  )
})

test_that("sve_from_model validates effect_name in data", {
  skip_if_not_installed("survival")

  data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    status = c(1, 1, 0, 1, 0, 1, 1, 0, 1, 0),
    treatment = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0),
    age = c(50, 55, 60, 65, 70, 52, 58, 62, 68, 72)
  )

  wrong_data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    status = c(1, 1, 0, 1, 0, 1, 1, 0, 1, 0),
    vaccine = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0),
    age = c(50, 55, 60, 65, 70, 52, 58, 62, 68, 72)
  )

  fit <- survival::coxph(survival::Surv(time, status) ~ treatment + age, data = data)

  expect_error(
    sve_from_model(fit, "treatment", data = wrong_data, method = "profile"),
    "not found in provided data"
  )
})

test_that("sve_from_model works with Cox model and profile method", {
  skip_if_not_installed("survival")

  set.seed(1)
  n <- 200
  data <- data.frame(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 60, 10)
  )

  fit <- survival::coxph(survival::Surv(time, status) ~ treatment + age, data = data)

  result <- sve_from_model(fit, "treatment", data = data, method = "profile")

  expect_true(is.data.frame(result))
  expect_equal(result$method, "Profile")
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$lower))
  expect_true(is.finite(result$upper))
  expect_lt(result$lower, result$upper)
})

test_that("sve_from_model works with Poisson GLM and profile method", {
  set.seed(1)
  n <- 200
  data <- data.frame(
    events = rpois(n, 5),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 60, 10)
  )

  fit <- glm(events ~ treatment + age, family = poisson(), data = data)

  result <- sve_from_model(fit, "treatment", data = data, method = "profile")

  expect_true(is.data.frame(result))
  expect_equal(result$method, "Profile")
  expect_true(is.finite(result$estimate))
  expect_lt(result$lower, result$upper)
})

test_that("sve_from_model Wald methods don't require data", {
  skip_if_not_installed("survival")

  set.seed(1)
  n <- 200
  data <- data.frame(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 60, 10)
  )

  fit <- survival::coxph(survival::Surv(time, status) ~ treatment + age, data = data)

  # Should work without data argument
  result_wald <- sve_from_model(fit, "treatment", method = "wald")
  result_tanh <- sve_from_model(fit, "treatment", method = "tanh-wald")

  expect_true(is.data.frame(result_wald))
  expect_true(is.data.frame(result_tanh))
  expect_equal(result_wald$estimate, result_tanh$estimate)
})

test_that("refit_model_constrained works for Cox models", {
  skip_if_not_installed("survival")

  set.seed(1)
  n <- 100
  data <- data.frame(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 60, 10)
  )

  fit <- survival::coxph(survival::Surv(time, status) ~ treatment + age, data = data)

  # Refit with treatment coefficient constrained to log(0.5)
  ll_constrained <- refit_model_constrained(model = fit,
                                            data = data,
                                            effect_name = "treatment",
                                            beta_fixed = log(0.5))
  ll_original <- as.numeric(logLik(fit))

  expect_true(is.finite(ll_constrained))
  expect_true(ll_constrained <= ll_original + 1e-6)
})

test_that("refit_model_constrained works for GLM models", {
  set.seed(1)
  n <- 100
  data <- data.frame(
    events = rpois(n, 5),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 60, 10)
  )

  fit <- glm(events ~ treatment + age, family = poisson(), data = data)

  # Refit with treatment coefficient constrained
  ll_constrained <- refit_model_constrained(fit, data, "treatment", log(0.8))
  ll_original <- as.numeric(logLik(fit))

  expect_true(is.finite(ll_constrained))
  expect_true(ll_constrained <= ll_original + 1e-6)
})
