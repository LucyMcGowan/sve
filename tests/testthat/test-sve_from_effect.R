test_that("sve_from_effect validates inputs", {
  expect_error(
    sve_from_effect(theta = -0.5, var_log_theta = 0.01),
    "must be positive"
  )
  expect_error(
    sve_from_effect(theta = 0.7, var_log_theta = -0.01),
    "must be non-negative"
  )
  expect_error(
    sve_from_effect(theta = c(0.7, 0.8), var_log_theta = 0.01),
    "must have the same length"
  )
})

test_that("sve_from_effect computes correct SVE transformations", {
  # Protective effect: theta < 1 -> SVE = 1 - theta
  expect_equal(
    sve_from_effect(theta = 0.5, var_log_theta = 0.04)$estimate,
    0.5
  )

  # Null effect: theta = 1 -> SVE = 0
  expect_equal(
    sve_from_effect(theta = 1.0, var_log_theta = 0.01)$estimate,
    0
  )

  # Harmful effect: theta > 1 -> SVE = (1 - theta) / theta
  expect_equal(
    sve_from_effect(theta = 1.5, var_log_theta = 0.04)$estimate,
    -1/3,
    tolerance = 1e-10
  )
})

test_that("sve_from_effect profile CI is properly ordered and bounded", {
  result <- sve_from_effect(
    theta = 0.7,
    var_log_theta = 0.0225,
    method = "profile"
  )

  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
  expect_gte(result$lower, -1)
  expect_lte(result$upper, 1)
})

test_that("sve_from_effect tanh-wald keeps bounds strictly in (-1, 1)", {
  # Strong effect with large variance - tanh transform ensures bounds
  result <- sve_from_effect(
    theta = 0.2,
    var_log_theta = 0.5,
    method = "tanh-wald"
  )

  expect_gt(result$lower, -1)
  expect_lt(result$upper, 1)
})


