#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption at steady state, with lag time
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V Central volume of distribution (L)
#' @param dur Duration of zero-order absorption (h)
#' @param dose Steady state dose
#' @param tau Dosing interval (h)
#' @param tlag Lag time (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) after dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_ss_1cmt_linear_oral_0_lag(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24, tlag=1.5)
#'
#' @export

calc_ss_1cmt_linear_oral_0_lag <-
  function(tad, CL, V, dur, dose, tau, tlag) {
    ### microconstants
    k   <- CL / V

    ### C(t) at steady state - eq 1.26 p. 14

    Ct <-
      (dose / dur) * (1 / (k * V)) * ((1 - exp(-k * dur)) * exp(-k * (tad - tlag - dur)) / (1 - exp(-k * tau)))

    Ct[tad < tlag] <- (dose / dur) * (1 / (k * V)) * ((1 - exp(-k * dur)) * exp(-k * (tad[tad < tlag] + tau - tlag - dur)) / (1 - exp(-k * tau)))

    Ct[tad >= tlag & tad < dur] <-
      (dose / dur) * (1 / (k * V)) * ((1 - exp(-k * (tad[tad >= tlag & tad < dur] - tlag))) +
                                        exp(-k * tau) * ((1 - exp(-k * dur)) * exp(-k * (tad[tad >= tlag & tad < dur] - tlag - dur)) /
                                                           (1 - exp(-k * tau))))

    Ct
  }
