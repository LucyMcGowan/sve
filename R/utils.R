check_proportions <- function(p0, p1) {
  if (any(p0 < 0 | p0 > 1)) {
    cli::cli_abort(c("x" =  "{.var p0} must be a proportion between 0 and 1."))
  }
  if (any(p1 < 0 | p1 > 1)) {
    cli::cli_abort(c("x" = "{.var p1} must be a proportion between 0 and 1."))
  }
  if (any(p0 == 0 & p1 == 0)) {
    cli::cli_abort(c("x" = "{.var p0} and {.var p1} are exactly 0.
                     Relative risk is not defined."))
  }
}

check_count_inputs <- function(x0, x1, n0, n1) {
  if (any(x0 < 0 | x1 < 0)) {
    cli::cli_abort(c("x" = "{.var x0} and {.var x1} must be non-negative."))
  }
  if (any(n0 <= 0 | n1 <= 0)) {
    cli::cli_abort(c("x" = "{.var n0} and {.var n1} must be positive."))
  }
  if (any(x0 > n0)) {
    cli::cli_abort(c("x" = "{.var x0} must not exceed {.var n0}."))
  }
  if (any(x1 > n1)) {
    cli::cli_abort(c("x" = "{.var x1} must not exceed {.var n1}."))
  }
}

check_confidence_level <- function(level) {
  if (any(level <= 0 | level >= 1)) {
    cli::cli_abort(c("x" = "{.var level} must be between 0 and 1."))
  }
}

check_theta <- function(theta, var_log_theta) {
  if (any(theta <= 0)) {
    cli::cli_abort("{.arg theta} must be positive.")
  }

  if (any(var_log_theta < 0)) {
    cli::cli_abort("{.arg var_log_theta} must be non-negative.")
  }

  if (length(theta) != length(var_log_theta)) {
    cli::cli_abort(
      c("{.arg theta} and {.arg var_log_theta} must have the same length.",
        "i" = "{.arg theta} has length {length(theta)}.",
        "i" = "{.arg var_log_theta} has length {length(var_log_theta)}.")
    )
  }
}
