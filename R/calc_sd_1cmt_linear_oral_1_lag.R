#' Calculate C(t) for a 1-compartment linear model with first-order absorption after a single oral dose, with lag time
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param ka First order absorption rate constant (/h)
#' @param tlag Lag time (h)
#' @param dose Dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_sd_1cmt_linear_oral_1_lag(t=0:24, CL=6, V=25, ka=1.1, dose=600, tlag=2)
#'
#' @export

calc_sd_1cmt_linear_oral_1_lag <- function(t, CL, V, ka, tlag, dose) {

  ### microconstants
  k   <- CL/V

  ### C(t) after single dose - eq 1.11 p. 9

  Ct <- (dose/V) * (ka / (ka - k)) * (exp(-k * (t - tlag)) - exp(-ka * (t - tlag)))
  Ct[t < tlag] <- 0

  Ct
}

