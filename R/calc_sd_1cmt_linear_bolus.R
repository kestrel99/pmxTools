#' Calculate C(t) for a 1-compartment linear model after a single IV bolus dose
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param dose Dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_sd_1cmt_linear_bolus(t=0:24, CL=6, V=25, dose=600)
#'
#' @export

calc_sd_1cmt_linear_bolus <- function(t, CL, V, dose) {

  ### microconstants
  k   <- CL/V

  ### C(t) after single dose - eq 1.1 p. 5

  Ct <- (dose/V) * exp(-k*t)

  Ct
}

