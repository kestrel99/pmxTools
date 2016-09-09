#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption after a single dose
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param dur Duration of zero-order absorption (h)
#' @param Dose Dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- sd1CptLinearOral0O(t=0:24, CL=6, V=25, dur=1.5, Dose=600)
#'

sd1CptLinearOral0O <- function(t, CL, V, dur, Dose) {

  ### microconstants 
  k   <- CL/V

  ### C(t) after single dose - eq 1.21 p. 13
  
  Ct <- (Dose / dur) * (1 / (k * V)) * (1 - exp(-k * dur)) * exp(-k * (t - dur))
  
  Ct[t <= dur] <- (Dose / dur) * (1 / (k * V)) * (1 - exp(-k * t[t <= dur]))

  Ct
}

