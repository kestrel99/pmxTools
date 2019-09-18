#' Calculate C(t) for a 1-compartment linear model at steady-state
#'
#' @param tad Time after dose (h)
#' @param tau Dosing interval (h)
#' @param dose Dose
#' @param dur Duration of zero-order absorption (h)
#' @param tinf Duration of infusion (h)
#' @param ... Passed to `calc_derived_1cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_ss_1cmt <- function(tad, tau, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with IV bolus dosing at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_bolus(t=0:24, CL=6, V=25, dose=600, tau=24)
#' @export
calc_ss_1cmt_linear_bolus <- function(tad, CL, V, dose, tau) {
  ### microconstants
  k   <- CL/V
  ### C(t) after single dose - eq 1.1 p. 5
  Ct <- (dose/V) * (exp(-k*tad))/(1 - exp(-k * tau))
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with infusion at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_infusion(tad=0:36, CL=2, V=25, dose=600, tinf=1, tau=24)
#' @export
calc_ss_1cmt_linear_infusion <- function(tad, CL, V, dose, tinf, tau) {
  ### microconstants
  k   <- CL/V
  ### C(t) after single dose - eq 1.8 p. 8
  Ct <- (dose/tinf) * (1 / (k*V) ) * ((1 - exp(-k * tinf)) * exp(-k * (tad - tinf))) / (1 - exp(-k * tau))
  Ct[tad <= tinf] <- (dose/tinf) * (1 / (k*V) ) * ((1 - exp(-k * tad[tad<=tinf])) +
                                                     exp(-k * tau) * (((1 - exp(-k * tinf)) * exp(-k * (tad[tad<=tinf] - tinf))) /
                                                                        (1 - exp(-k * tau))))
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_0(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24)
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

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption at steady state, with lag time
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_0_lag(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24, tlag=1.5)
#' @export
calc_ss_1cmt_linear_oral_0_lag <- function(tad, CL, V, dur, dose, tau, tlag) {
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

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state, with lag time
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_1_lag(tad=0:36, CL=2, V=25, dose=600,
#'     ka=0.25, tlag=0.75, tau=24)
#' @export
calc_ss_1cmt_linear_oral_1_lag <- function(tad, CL, V, ka, dose, tlag, tau) {
  ### microconstants
  k   <- CL / V
  ### C(t) after single dose - eq 1.16 p. 11
  Ct <-
    (dose / V) * (ka / (ka - k)) * ((exp(-k * (tad - tlag)) / (1 - exp(-k * tau))) - (exp(-ka * (tad - tlag)) / (1 - exp(-ka * tau))))
  
  Ct[tad < tlag] <-
    (dose / V) * (ka / (ka - k)) * ((exp(-k * (tad[tad < tlag] + tau - tlag)) / (1 - exp(-k * tau))) - (exp(-ka * (tad[tad < tlag] + tau - tlag)) / (1 - exp(-ka * tau))))
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_1(tad=0:36, CL=2, V=25, dose=600, ka=0.25, tau=24)
#' @export
calc_ss_1cmt_linear_oral_1 <- function(tad, CL, V, ka, dose, tau) {
  ### microconstants
  k   <- CL/V
  ### C(t) after single dose - eq 1.8 p. 8
  Ct <- (dose/V) * (ka / (ka - k)) * ((exp(-k * tad) / (1 - exp(-k * tau))) - (exp(-ka * tad) / (1 - exp(-ka * tau))))
  Ct
}
