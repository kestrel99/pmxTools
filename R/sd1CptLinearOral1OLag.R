#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption after a single dose, with lag time
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param ka First order absorption rate constant (/h)
#' @param tlag Lag time (h)
#' @param Dose Dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- sd1CptLinearOral1O(t=0:24, CL=6, V=25, ka=1.1, Dose=600, tlag=2)
#'

sd1CptLinearOral1OLag <- function(t, CL, V, ka, tlag, Dose) {

  ### microconstants 
  k   <- CL/V

  ### C(t) after single dose - eq 1.11 p. 9

  Ct <- (Dose/V) * (ka / (ka - k)) * (exp(-k * (t - tlag)) - exp(-ka * (t - tlag)))
  Ct[t < tlag] <- 0

  Ct
}

