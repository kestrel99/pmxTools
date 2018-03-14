#' Calculate C(t) for a 1-compartment linear model with IV bolus dosing at steady state
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param dose Steady state dose
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) after dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_ss_1cmt_linear_bolus(t=0:24, CL=6, V=25, dose=600, tau=24)
#'
#' @export

calc_ss_1cmt_linear_bolus <- function(tad, CL, V, dose, tau) {

  ### microconstants
  k   <- CL/V

  ### C(t) after single dose - eq 1.1 p. 5

  Ct <- (dose/V) * (exp(-k*tad))/(1 - exp(-k * tau))

  Ct
}

