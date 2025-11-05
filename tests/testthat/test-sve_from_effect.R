test_that("sve_from_effect validates inputs", {
  expect_error(sve_from_effect(theta = -0.5, var_log_theta = 0.01))
  expect_error(sve_from_effect(theta = 0.7, var_log_theta = -0.01))
  expect_error(sve_from_effect(theta = c(0.7, 0.8), var_log_theta = 0.01))
})

test_that("sve_from_effect calculates correct point estimates", {
  # theta = 0.5: SVE = 0.5
  result <- sve_from_effect(theta = 0.5, var_log_theta = 0.04)
  expect_equal(result$estimate, 0.5)

  # theta = 1: SVE = 0
  result <- sve_from_effect(theta = 1.0, var_log_theta = 0.01)
  expect_equal(result$estimate, 0)

  # theta = 1.5: SVE = -1/3
  result <- sve_from_effect(theta = 1.5, var_log_theta = 0.04)
  expect_equal(result$estimate, -1/3, tolerance = 1e-10)
})

test_that("sve_from_effect handles vectorized inputs", {
  hrs <- c(0.5, 0.8, 1.2)
  var_log_hrs <- c(0.04, 0.025, 0.036)

  result <- sve_from_effect(theta = hrs, var_log_theta = var_log_hrs)

  expect_equal(nrow(result), 3)
  expect_equal(result$estimate[1], 0.5)
  expect_gt(result$estimate[2], 0)
  expect_lt(result$estimate[3], 0)
})

test_that("sve_from_effect confidence intervals are properly ordered", {
  result <- sve_from_effect(theta = 0.7, var_log_theta = 0.0225)
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("sve_from_effect tanh-wald keeps bounds in valid range", {
  # Strong effect with large variance
  result <- sve_from_effect(theta = 0.2, var_log_theta = 0.5,
                            method = "tanh-wald")
  expect_gte(result$lower, -1)
  expect_lte(result$upper, 1)
})
