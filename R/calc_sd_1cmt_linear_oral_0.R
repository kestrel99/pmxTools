#' Calculate C(t) for a 1-compartment linear model with zero-order absorption after a single oral dose
#'
#' @param t Time after dose (h)
#' @param ... Passed to `calc_derived_1cpt()`
#' @param dur Duration of zero-order absorption (h)
#' @param dose Dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#'
#' @examples
#' Ctrough <- calc_sd_1cmt_linear_oral_0(t=0:24, CL=6, V=25, dur=1.5, dose=600)
#' @export
calc_sd_1cmt_linear_oral_0 <- function(t, ..., dur, dose) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.21 p. 13
  Ct <- (dose/ dur) * (1 / (param$k10 * param$V1)) * (1 - exp(-param$k10 * dur)) * exp(-param$k10 * (t - dur))
  Ct[t <= dur ] <- (dose/ dur) * (1 / (param$k10 * param$V1)) * (1 - exp(-param$k10 * t[t <= dur]))
  Ct
}
