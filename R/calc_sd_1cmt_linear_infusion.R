#' Calculate C(t) for a 1-compartment linear model after a single IV infusion
#'
#' @param t Time after dose (h)
#' @param ... Passed to `calc_derived_1cpt()`
#' @param dose Dose
#' @param tinf Duration of infusion (h)
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#'
#' @examples
#' Ctrough <- calc_sd_1cmt_linear_infusion(t=0:24, CL=6, V=25, dose=600, tinf=1)
#' @export
calc_sd_1cmt_linear_infusion <- function(t, ..., dose, tinf) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.6 p. 7
  Ct <- (dose/tinf) * (1 / (param$k10*param$V1) ) * (1 - exp(-param$k10 * tinf)) * exp(-param$k10 * (t-tinf))
  Ct[t <= tinf] <- (dose/tinf) * (1 / (param$k10*param$V1) ) * (1 - exp(-param$k10 * t[t <= tinf]))
  Ct
}
