#' Calculate C(t) for a 1-compartment linear model with IV bolus at steady state
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param Dose Steady state dose
#' @param Tinf Duration of infusion (h)
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) after dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- ss1CptLinearInf(tad=0:36, CL=2, V=25, Dose=600, Tinf=1, tau=24)
#'

ss1CptLinearInf <- function(tad, CL, V, Dose, Tinf, tau) {

  ### microconstants 
  k   <- CL/V

  ### C(t) after single dose - eq 1.8 p. 8
  
  Ct <- (Dose/Tinf) * (1 / (k*V) ) * ((1 - exp(-k * Tinf)) * exp(-k * (tad - Tinf))) / (1 - exp(-k * tau))
  
  Ct[tad <= Tinf] <- (Dose/Tinf) * (1 / (k*V) ) * ((1 - exp(-k * tad[tad<=Tinf])) +  
                                                     exp(-k * tau) * (((1 - exp(-k * Tinf)) * exp(-k * (tad[tad<=Tinf] - Tinf))) / 
                                                     (1 - exp(-k * tau))))

  Ct
}

