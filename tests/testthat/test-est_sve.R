test_that("est_sve validates inputs", {
  expect_error(
    est_sve(x0 = -1, x1 = 50, n0 = 1000, n1 = 1000),
    "non-negative"
  )
  expect_error(
    est_sve(x0 = 100, x1 = 50, n0 = 50, n1 = 1000),
    "x0.*"
  )
})

test_that("est_sve computes correct SVE transformations", {
  # Protective effect: p0 > p1
  result <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000)
  expect_equal(result$estimate, 0.5)

  # Null effect: p0 = p1
  result_null <- est_sve(x0 = 100, x1 = 100, n0 = 1000, n1 = 1000)
  expect_equal(result_null$estimate, 0)

  # Harmful effect: p0 < p1
  result_harm <- est_sve(x0 = 50, x1 = 100, n0 = 1000, n1 = 1000)
  expect_equal(result_harm$estimate, -0.5)
})


test_that("est_sve all methods produce same point estimate", {
  x0 <- 100; x1 <- 50; n0 <- 1000; n1 <- 1000

  estimates <- c(
    est_sve(x0, x1, n0, n1, method = "profile")$estimate,
    est_sve(x0, x1, n0, n1, method = "tanh-wald")$estimate,
    est_sve(x0, x1, n0, n1, method = "wald")$estimate,
    est_sve(x0, x1, n0, n1, method = "exact")$estimate
  )

  expect_true(all(abs(diff(estimates)) < 1e-10))
})

test_that("est_sve profile method produces valid confidence intervals", {
  result <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000, method = "profile")

  expect_equal(result$method, "Profile")
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
  expect_gt(result$lower, -1)
  expect_lt(result$upper, 1)
})

test_that("est_sve profile method handles small samples", {
  result <- est_sve(x0 = 10, x1 = 5, n0 = 100, n1 = 100, method = "profile")

  expect_true(is.finite(result$lower))
  expect_true(is.finite(result$upper))
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("est_sve tanh-wald keeps bounds strictly in (-1, 1)", {
  # Extreme case with small sample
  result <- est_sve(x0 = 10, x1 = 1, n0 = 20, n1 = 20, method = "tanh-wald")

  expect_gt(result$lower, -1)
  expect_lt(result$upper, 1)
})

test_that("est_sve smooth parameter controls variance smoothing", {
  # Near null case where smoothing matters
  result_smooth <- est_sve(x0 = 100, x1 = 95, n0 = 1000, n1 = 1000,
                           method = "wald", smooth = TRUE)
  result_no_smooth <- est_sve(x0 = 100, x1 = 95, n0 = 1000, n1 = 1000,
                              method = "wald", smooth = FALSE)

  # Smoothing should affect the CI width near null
  expect_false(isTRUE(all.equal(result_smooth$lower, result_no_smooth$lower)))
})

test_that("est_sve warns when epsilon provided but smooth = FALSE", {
  expect_warning(
    est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000,
            method = "wald", smooth = FALSE, epsilon = 0.1),
    "epsilon.*ignored.*smooth = FALSE"
  )
})

