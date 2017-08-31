#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param ka First order absorption rate constant (/h)
#' @param dose Steady state dose
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) after dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_ss_1cmt_linear_oral_1(tad=0:36, CL=2, V=25, dose=600, ka=0.25, tau=24)
#'
#' @export

calc_ss_1cmt_linear_oral_1 <- function(tad, CL, V, ka, dose, tau) {

  ### microconstants
  k   <- CL/V

  ### C(t) after single dose - eq 1.8 p. 8

  Ct <- (dose/V) * (ka / (ka - k)) * ((exp(-k * tad) / (1 - exp(-k * tau))) - (exp(-ka * tad) / (1 - exp(-ka * tau))))

  Ct
}

