#' Calculate C(t) for a 2-compartment linear model after a single zero-order oral dose
#'
#' @param t Time after dose (h)
#' @param ... Passed to `calc_derived_2cpt()`
#' @param dur Duration of zero-order absorption (h)
#' @param dose Steady state dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_sd_2cmt_linear_oral_0(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1)
#' @export
calc_sd_2cmt_linear_oral_0 <- function(t, ..., dur, dose) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.4 p. 28
  A <- (1/param$V1) * ((param$alpha - param$k21) / (param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21) / (param$beta - param$alpha))
  ### C(t) after single dose - eq 1.51 p. 29
  Ct <- (dose / dur) * (((A / param$alpha) * (1 - exp(-param$alpha * dur)) * exp(-param$alpha * (t - dur))) +
                      ((B / param$beta) * (1 - exp(-param$beta * dur)) * exp(-param$beta * (t - dur))))
  Ct[t < dur] <- (dose / dur) * ((A / param$alpha) * (1 - exp(-param$alpha * t[t < dur])) +
                                 (B / param$beta) * (1 - exp(-param$beta * t[t < dur])))
  Ct
}
