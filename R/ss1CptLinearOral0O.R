#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param dur Duration of zero-order absorption (h)
#' @param Dose Steady state dose
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) after dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- ss1CptLinearOral1O(tad=0:36, CL=2, V=25, Dose=600, ka=0.25, tau=24)
#'

ss1CptLinearOral0O <- function(tad, CL, V, dur, Dose, tau) {

  ### microconstants 
  k   <- CL/V

  ### C(t) at steady state - eq 1.23 p. 14
  
  Ct <- (Dose / dur) * (1 / (k * V)) * ((1- exp(-k * dur)) * exp(-k * (tad - dur)) / (1 - exp(-k * tau)))
  
  Ct[tad <= dur] <- (Dose / dur) * (1 / (k * V)) * ( (1 - exp(-k * tad[tad <= dur])) + 
                                                       exp(-k * tau)*( (1 - exp(-k * dur)) * exp(-k * (tad[tad <= dur] - dur)) / 
                                                                         (1 - exp(-k * tau))))

  Ct
}

