#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state, with lag time
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param ka First order absorption rate constant (/h)
#' @param dose Steady state dose
#' @param tlag Lag time (h)
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
#' Ctrough <- calc_ss_1cmt_linear_oral_1_lag(tad=0:36, CL=2, V=25, dose=600,
#'     ka=0.25, tlag=0.75, tau=24)
#'
#' @export

calc_ss_1cmt_linear_oral_1_lag <-
  function(tad, CL, V, ka, dose, tlag, tau) {
    ### microconstants
    k   <- CL / V

    ### C(t) after single dose - eq 1.16 p. 11

    Ct <-
      (dose / V) * (ka / (ka - k)) * ((exp(-k * (tad - tlag)) / (1 - exp(-k * tau))) - (exp(-ka * (tad - tlag)) / (1 - exp(-ka * tau))))

    Ct[tad < tlag] <-
      (dose / V) * (ka / (ka - k)) * ((exp(-k * (tad[tad < tlag] + tau - tlag)) / (1 - exp(-k * tau))) - (exp(-ka * (tad[tad < tlag] + tau - tlag)) / (1 - exp(-ka * tau))))

    Ct
  }
