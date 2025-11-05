test_that("est_sve calculates correct point estimates", {
  # When p0 = 0.1 and p1 = 0.05, SVE = 0.5
  result <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000)
  expect_equal(result$estimate, 0.5)

  # When proportions are equal, SVE should be 0
  result_null <- est_sve(x0 = 100, x1 = 100, n0 = 1000, n1 = 1000)
  expect_equal(result_null$estimate, 0)

  # Harmful effect (negative SVE)
  result_harm <- est_sve(x0 = 50, x1 = 100, n0 = 1000, n1 = 1000)
  expect_equal(result_harm$estimate, -0.5)
})

test_that("est_sve confidence intervals are properly ordered", {
  result <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000)
  expect_lt(result$lower, result$estimate)
  expect_gt(result$upper, result$estimate)
})

test_that("est_sve handles boundary cases", {
  # Perfect efficacy
  result <- est_sve(x0 = 100, x1 = 0, n0 = 1000, n1 = 1000)
  expect_equal(result$estimate, 1)

  # Perfect harm
  result <- est_sve(x0 = 0, x1 = 100, n0 = 1000, n1 = 1000)
  expect_equal(result$estimate, -1)

  # Exact method with zero events
  result <- est_sve(x0 = 50, x1 = 0, n0 = 100, n1 = 100, method = "exact")
  expect_true(is.finite(result$lower))
  expect_true(is.finite(result$upper))
})

test_that("est_sve all methods give same point estimate", {
  result_tanh <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000,
                         method = "tanh-wald")
  result_wald <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000,
                         method = "wald")
  result_exact <- est_sve(x0 = 100, x1 = 50, n0 = 1000, n1 = 1000,
                          method = "exact")

  expect_equal(result_tanh$estimate, result_wald$estimate)
  expect_equal(result_tanh$estimate, result_exact$estimate)
})
