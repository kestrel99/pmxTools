#' Calculate C(t) for a 2-compartment linear model after a single IV bolus dose
#'
#' @param t Time after dose (h)
#' @param ... Passed to `calc_derived_2cpt()`
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
#' Ctrough <- calc_sd_2cmt_linear_bolus(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10)
#' @export
calc_sd_2cmt_linear_bolus <- function(t, ..., dose) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.2 p. 21
  A <- (1/param$V1) * ((param$alpha - param$k21)/(param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21)/(param$beta - param$alpha))
  ### C(t) after single dose - eq 1.36 p. 22
  Ct <- dose * (A * exp(-param$alpha * t) + (B * exp(-param$beta * t)))
  Ct
}
