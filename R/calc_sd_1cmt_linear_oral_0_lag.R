#' Calculate C(t) for a 1-compartment linear model with zero-order absorption after a single oral dose, with lag time
#'
#' @param t Time after dose (h)
#' @inheritParams calc_derived
#' @param dur Duration of zero-order absorption (h)
#' @param dose Dose
#' @param tlag Lag time (h)
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_sd_1cmt_linear_oral_0_lag(t=0:24, CL=6, V=25, dur=1.5, dose=600, tlag=1.5)
#'
#' @export

calc_sd_1cmt_linear_oral_0_lag <- function(t, CL, V, dur, dose, tlag) {

  ### microconstants
  k   <- CL/V

  ### C(t) after single dose - eq 1.24 p. 13

  Ct <- (dose / dur) * (1 / (k * V)) * (1 - exp(-k * dur)) * exp(-k * (t - tlag - dur))

  Ct[t < tlag] <- 0

  Ct[t >= tlag & t < dur] <- (dose / dur) * (1 / (k * V)) * (1 - exp(-k * (t[t >= tlag & t < dur] - tlag)))

  Ct
}

