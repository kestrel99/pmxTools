#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption at steady state
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param dur Duration of zero-order absorption (h)
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
#' Ctrough <- calc_ss_1cmt_linear_oral_0(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24)
#'
#' @export

calc_ss_1cmt_linear_oral_0 <- function(tad, CL, V, dur, dose, tau) {

  ### microconstants
  k   <- CL/V

  ### C(t) at steady state - eq 1.23 p. 14

  Ct <- (dose / dur) * (1 / (k * V)) * ((1- exp(-k * dur)) * exp(-k * (tad - dur)) / (1 - exp(-k * tau)))

  Ct[tad <= dur] <- (dose / dur) * (1 / (k * V)) * ( (1 - exp(-k * tad[tad <= dur])) +
                                                       exp(-k * tau)*( (1 - exp(-k * dur)) * exp(-k * (tad[tad <= dur] - dur)) /
                                                                         (1 - exp(-k * tau))))

  Ct
}

