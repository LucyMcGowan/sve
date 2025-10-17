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
