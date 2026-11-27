# Choose the dosing regimen implied by the supplied absorption/infusion
# parameters, matching the suffixes of the calc_sd_*()/calc_ss_*() functions.
regimen_name <- function(dur = NULL, tinf = NULL, ka = NULL, tlag = NULL) {
  if (!is.null(tinf)) {
    if (!is.null(dur) || !is.null(ka) || !is.null(tlag)) {
      stop("`tinf` (IV infusion) cannot be combined with `dur`, `ka` or `tlag`.", call. = FALSE)
    }
    return("infusion")
  }
  if (!is.null(dur) && !is.null(ka)) {
    stop("Specify only one of `dur` (zero-order absorption) or `ka` (first-order absorption).", call. = FALSE)
  }
  absorption <- if (!is.null(dur)) "oral_0" else if (!is.null(ka)) "oral_1"
  if (is.null(absorption)) {
    if (!is.null(tlag)) stop("`tlag` requires oral absorption (`ka` or `dur`).", call. = FALSE)
    return("bolus")
  }
  if (is.null(tlag)) absorption else paste0(absorption, "_lag")
}

# Call the regimen-specific function `<prefix>_linear_<regimen>()`.
dispatch_regimen <- function(prefix, args, dur = NULL, tinf = NULL, ...) {
  dots <- list(...)
  regimen <- regimen_name(dur = dur, tinf = tinf, ka = dots[["ka"]], tlag = dots[["tlag"]])
  fun <- get(paste0(prefix, "_linear_", regimen), mode = "function")
  do.call(fun, c(args, Filter(Negate(is.null), list(dur = dur, tinf = tinf)), dots))
}
