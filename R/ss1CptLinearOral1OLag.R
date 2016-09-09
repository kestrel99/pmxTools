#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state, with lag time
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param ka First order absorption rate constant (/h)
#' @param Dose Steady state dose
#' @param tlag Lag time (h)
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) after dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- ss1CptLinearOral1O(tad=0:36, CL=2, V=25, Dose=600, ka=0.25, tlag=0.75, tau=24)
#'

ss1CptLinearOral1OLag <- function(tad, CL, V, ka, Dose, tlag, tau) {

  ### microconstants 
  k   <- CL/V

  ### C(t) after single dose - eq 1.16 p. 11
  
  Ct <- (Dose/V) * (ka / (ka - k)) * ((exp(-k * (tad - tlag)) / (1 - exp(-k * tau))) - (exp(-ka * (tad - tlag)) / (1 - exp(-ka * tau))))
  
  Ct[tad < tlag] <- (Dose/V) * (ka / (ka - k)) * ((exp(-k * (tad[tad < tlag] + tau - tlag)) / (1 - exp(-k * tau))) - (exp(-ka * (tad[tad < tlag] + tau - tlag)) / (1 - exp(-ka * tau))))
  
  Ct
}

