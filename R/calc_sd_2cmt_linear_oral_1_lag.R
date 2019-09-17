#' Calculate C(t) for a 2-compartment linear model after a single first-order oral dose with a lag time
#'
#' @param t Time after dose (h)
#' @param ... Passed to `calc_derived_2cpt()`
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
#' Ctrough <- calc_sd_2cmt_linear_oral_1_lag(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tlag = 2)
#' @export
calc_sd_2cmt_linear_oral_1_lag <- function(t, ..., dose) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.3 p. 24
  A <- (param$ka/param$V1) * ((param$k21 - param$alpha) / ((param$ka - param$alpha) * (param$beta - param$alpha)))
  B <- (param$ka/param$V1) * ((param$k21 - param$beta) / ((param$ka - param$beta) * (param$alpha - param$beta)))
  ### C(t) after single dose - eq 1.41 p. 25
  Ct <- dose * ((A * exp(-param$alpha * (t - param$tlag))) + (B * exp(-param$beta * (t - param$tlag))) - ((A + B) * exp(-param$ka * (t - param$tlag))))
  Ct[t < param$tlag] <- 0
  Ct
}
