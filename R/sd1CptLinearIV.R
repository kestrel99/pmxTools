#' Calculate C(t) for a 1-compartment linear model with IV bolus after a single dose
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param Dose Dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- sd1CptLinearIV(t=0:24, CL=6, V=25, Dose=600)
#'

sd1CptLinearIV <- function(t, CL, V, Dose) {

  ### microconstants 
  k   <- CL/V

  ### C(t) after single dose - eq 1.1 p. 5

  Ct <- (Dose/V) * exp(-k*t)

  Ct
}

