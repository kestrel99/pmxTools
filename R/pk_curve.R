#' Provide concentration-time curves.
#'
#' Concentrations after repeated dosing are obtained by superposition: each dose
#' given at time \eqn{n \cdot ii} (\eqn{n = 0, \ldots, addl}) contributes the
#' single-dose concentration at the time elapsed since that dose, and doses not
#' yet given contribute nothing. Observation times need not be regularly spaced
#' or start at zero.
#'
#' @param t Observation time in h, specified as a vector of non-negative values
#'   relative to the first dose.
#' @param model The model to use. Must be one of "1cmt_bolus", "1cmt_infusion", "1cmt_oral", "2cmt_bolus", "2cmt_infusion",
#' "2cmt_oral", "3cmt_bolus", "3cmt_infusion", "3cmt_oral". The default is "1cmt_oral".
#' @param params A named list containing parameter values for the selected model type.
#' @param dose Dose amount.
#' @param ii Interdose interval (or tau), in hours (default 24).
#' @param addl Number of additional doses (default 0).
#' @param ss Assume steady state concentration (default \code{FALSE}). When
#'   \code{TRUE}, steady state is assumed to have been reached with dosing
#'   interval \code{ii} before the dose at time 0 (as for \code{SS=1} in
#'   NONMEM), and the \code{addl} additional doses are given on top of it. After
#'   the last dose the curve declines as the drug washes out.
#'
#' @return A data frame containing times (\code{t}) and concentrations (\code{cp}).
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @export
#'
#' @examples
#' plot(pk_curve(t=seq(0,72,by=0.1), model="3cmt_oral", ii=12, addl=5,
#'   params=list(CL=2.5, V1=25, V2=2, V3=5, Q2=0.5, Q3=0.25, ka=1)), type="l")
#'
#' # Steady state, followed by washout
#' plot(pk_curve(t=seq(0,72,by=0.1), model="1cmt_oral", ii=12, ss=TRUE,
#'   params=list(CL=2.5, V=25, ka=1)), type="l")
pk_curve <- function(t, model = "1cmt_oral", params=list(ka=2.77, CL=2.5, V=25), dose=600, ii=24, addl=0, ss=FALSE) {
  if (!is.character(model) || length(model) != 1 || !(model %in% names(pk_curve_models))) {
    stop("`model` must be one of: ", paste0('"', names(pk_curve_models), '"', collapse = ", "), call. = FALSE)
  }
  spec <- pk_curve_models[[model]]
  if (!setequal(names(params), spec$params)) {
    stop("The ", spec$label, " model requires parameter values for ",
         paste(spec$params, collapse = ", "), ".", call. = FALSE)
  }
  # Every dose must use the same parameter set, so each parameter is a single
  # value (the curve describes one individual)
  non_scalar <- names(params)[lengths(params) != 1]
  if (length(non_scalar) > 0) {
    stop("Each parameter in `params` must be a single value; check: ", paste(non_scalar, collapse = ", "),
         ". Call pk_curve() once per parameter set.", call. = FALSE)
  }
  assert_nonnegative(t)
  assert_positive(dose)
  assert_positive(ii)
  assert_nonnegative(addl)
  if (length(dose) != 1 || length(ii) != 1) {
    stop("`dose` and `ii` must be single values.", call. = FALSE)
  }
  if (length(addl) != 1 || addl != round(addl)) {
    stop("`addl` must be a single non-negative whole number.", call. = FALSE)
  }
  if (!is.logical(ss) || length(ss) != 1 || is.na(ss)) {
    stop("`ss` must be TRUE or FALSE.", call. = FALSE)
  }

  single_dose <- function(elapsed) do.call(spec$sd, c(list(t = elapsed, dose = dose), params))

  cp <- numeric(length(t))
  if (ss) {
    # Contribution of the dose at time 0 and every dose before it
    cp <- do.call(spec$ss, c(list(tad = t, tau = ii, dose = dose), params))
  }
  dose_times <- ii * seq(from = if (ss) 1 else 0, length.out = addl + !ss)
  for (dose_time in dose_times) {
    # Allow for rounding in dose_time (3 * 0.1 > 0.3), so a dose scheduled at an
    # observation time is counted there, at zero elapsed time
    elapsed <- t - dose_time
    given <- elapsed >= -sqrt(.Machine$double.eps) * dose_time
    if (any(given)) {
      cp[given] <- cp[given] + single_dose(pmax(elapsed[given], 0))
    }
  }

  data.frame(t = t, cp = cp)
}

pk_curve_models <- list(
  "1cmt_oral" = list(label = "one-compartment oral", params = c("CL", "V", "ka"),
                     sd = calc_sd_1cmt_linear_oral_1, ss = calc_ss_1cmt_linear_oral_1),
  "1cmt_bolus" = list(label = "one-compartment IV bolus", params = c("CL", "V"),
                      sd = calc_sd_1cmt_linear_bolus, ss = calc_ss_1cmt_linear_bolus),
  "1cmt_infusion" = list(label = "one-compartment IV infusion", params = c("CL", "V", "tinf"),
                         sd = calc_sd_1cmt_linear_infusion, ss = calc_ss_1cmt_linear_infusion),
  "2cmt_oral" = list(label = "two-compartment oral", params = c("CL", "V1", "V2", "Q", "ka"),
                     sd = calc_sd_2cmt_linear_oral_1, ss = calc_ss_2cmt_linear_oral_1),
  "2cmt_bolus" = list(label = "two-compartment IV bolus", params = c("CL", "V1", "V2", "Q"),
                      sd = calc_sd_2cmt_linear_bolus, ss = calc_ss_2cmt_linear_bolus),
  "2cmt_infusion" = list(label = "two-compartment IV infusion", params = c("CL", "V1", "V2", "Q", "tinf"),
                         sd = calc_sd_2cmt_linear_infusion, ss = calc_ss_2cmt_linear_infusion),
  "3cmt_oral" = list(label = "three-compartment oral", params = c("CL", "V1", "V2", "V3", "Q2", "Q3", "ka"),
                     sd = calc_sd_3cmt_linear_oral_1, ss = calc_ss_3cmt_linear_oral_1),
  "3cmt_bolus" = list(label = "three-compartment IV bolus", params = c("CL", "V1", "V2", "V3", "Q2", "Q3"),
                      sd = calc_sd_3cmt_linear_bolus, ss = calc_ss_3cmt_linear_bolus),
  "3cmt_infusion" = list(label = "three-compartment IV infusion", params = c("CL", "V1", "V2", "V3", "Q2", "Q3", "tinf"),
                         sd = calc_sd_3cmt_linear_infusion, ss = calc_ss_3cmt_linear_infusion)
)
