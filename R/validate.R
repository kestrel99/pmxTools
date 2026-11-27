# Input validation helpers shared by the PK calculation functions.
#
# Each helper is a no-op for NULL so optional arguments can be checked
# unconditionally, and stops with a message naming the offending argument.

check_numeric_domain <- function(x, name, lower, inclusive) {
  if (is.null(x)) return(invisible(x))
  ok_type <- is.numeric(x) && length(x) > 0 && !anyNA(x)
  ok_range <- ok_type && all(if (inclusive) x >= lower else x > lower)
  if (!ok_range) {
    bound <- if (inclusive) "non-negative" else "positive"
    stop("`", name, "` must be ", bound, " numeric value(s) with no missing values.", call. = FALSE)
  }
  invisible(x)
}

assert_positive <- function(x, name = deparse(substitute(x))) {
  check_numeric_domain(x, name, lower = 0, inclusive = FALSE)
}

assert_nonnegative <- function(x, name = deparse(substitute(x))) {
  check_numeric_domain(x, name, lower = 0, inclusive = TRUE)
}

# Validate the arguments that the calc_sd_*() and calc_ss_*() functions take
# directly (structural parameters are checked by calc_derived_*()).
check_conc_inputs <- function(dose, t = NULL, tad = NULL, tau = NULL, tinf = NULL, dur = NULL) {
  assert_nonnegative(t)
  assert_nonnegative(tad)
  assert_positive(tau)
  assert_nonnegative(dose)
  assert_positive(tinf)
  assert_positive(dur)
  if (!is.null(tau)) {
    if (!is.null(tinf) && any(tinf > tau)) {
      stop("`tinf` must not exceed `tau`: the steady-state solution does not cover overlapping infusions.", call. = FALSE)
    }
    if (!is.null(dur) && any(dur > tau)) {
      stop("`dur` must not exceed `tau`: the steady-state solution does not cover overlapping absorption.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# First-order absorption solutions divide by (ka - rate) for each disposition
# rate constant, so they are singular when ka coincides with one of them.
check_ka_distinct <- function(param) {
  tol <- sqrt(.Machine$double.eps)
  rates <- intersect(c("alpha", "beta", "gamma"), names(param))
  coincident <- vapply(rates, function(r) any(abs(param$ka - param[[r]]) <= tol * param[[r]]), logical(1))
  if (any(coincident)) {
    stop("`ka` equals the disposition rate constant ", paste(rates[coincident], collapse = ", "),
         ", where the closed-form solution is singular. Adjust `ka` slightly.", call. = FALSE)
  }
  invisible(TRUE)
}

# Validate named lists of optional (possibly NULL) inputs to calc_derived_*().
check_derived_inputs <- function(positive = list(), nonnegative = list()) {
  for (nm in names(positive)) assert_positive(positive[[nm]], nm)
  for (nm in names(nonnegative)) assert_nonnegative(nonnegative[[nm]], nm)
  invisible(TRUE)
}
