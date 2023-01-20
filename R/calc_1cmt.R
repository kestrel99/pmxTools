# Single-dose #####

#' Calculate C(t) for a 1-compartment linear model
#'
#' @param t Time after dose (h)
#' @param dose Dose
#' @param dur Duration of zero-order absorption (h)
#' @param tinf Duration of infusion (h)
#' @param ... Passed to `calc_derived_1cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @author Bill Denney, \email{wdenney@@humanpredictions.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_sd_1cmt <- function(t, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_sd_1cmt
#' Calculate C(t) for a 1-compartment linear model after a single IV bolus dose
#' @examples
#' Ct <- calc_sd_1cmt_linear_bolus(t=0:24, CL=6, V=25, dose=600)
#' @export
calc_sd_1cmt_linear_bolus <- function(t, dose, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.1 p. 5
  Ct <- (dose/param$V1) * exp(-param$k10*t)
  Ct
}

#' @describeIn calc_sd_1cmt
#' Calculate C(t) for a 1-compartment linear model with first-order absorption after a single oral dose, with lag time
#' @examples
#' Ct <- calc_sd_1cmt_linear_oral_1_lag(t=0:24, CL=6, V=25, ka=1.1, dose=600, tlag=2)
#' @export
calc_sd_1cmt_linear_oral_1_lag <- function(t, dose, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.11 p. 9
  Ct <- (dose/param$V1) * (param$ka / (param$ka - param$k10)) * (exp(-param$k10 * (t - param$tlag)) - exp(-param$ka * (t - param$tlag)))
  Ct[t < param$tlag] <- 0
  Ct
}

#' @describeIn calc_sd_1cmt
#' Calculate C(t) for a 1-compartment linear model after a single IV infusion
#' @examples
#' Ct <- calc_sd_1cmt_linear_infusion(t=0:24, CL=6, V=25, dose=600, tinf=1)
#' @export
calc_sd_1cmt_linear_infusion <- function(t, dose, tinf, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.6 p. 7
  Ct <- (dose/tinf) * (1 / (param$k10*param$V1) ) * (1 - exp(-param$k10 * tinf)) * exp(-param$k10 * (t-tinf))
  Ct[t <= tinf] <- (dose/tinf) * (1 / (param$k10*param$V1) ) * (1 - exp(-param$k10 * t[t <= tinf]))
  Ct
}

#' @describeIn calc_sd_1cmt
#' Calculate C(t) for a 1-compartment linear model with zero-order absorption after a single oral dose
#' @examples
#' Ct <- calc_sd_1cmt_linear_oral_0(t=0:24, CL=6, V=25, dur=1.5, dose=600)
#' @export
calc_sd_1cmt_linear_oral_0 <- function(t, dose, dur, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.21 p. 13
  Ct <- (dose/ dur) * (1 / (param$k10 * param$V1)) * (1 - exp(-param$k10 * dur)) * exp(-param$k10 * (t - dur))
  Ct[t <= dur ] <- (dose/ dur) * (1 / (param$k10 * param$V1)) * (1 - exp(-param$k10 * t[t <= dur]))
  Ct
}

#' @describeIn calc_sd_1cmt
#' Calculate C(t) for a 1-compartment linear model with first-order absorption after a single oral dose
#' @examples
#' Ct <- calc_sd_1cmt_linear_oral_1(t=0:24, CL=6, V=25, ka=1.1, dose=600)
#' @export
calc_sd_1cmt_linear_oral_1 <- function(t, dose, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.11 p. 9
  Ct <- (dose/param$V1) * (param$ka / (param$ka - param$k10)) * (exp(-param$k10 * t) - exp(-param$ka * t))
  Ct
}

#' @describeIn calc_sd_1cmt
#' Calculate C(t) for a 1-compartment linear model with zero-order absorption after a single oral dose, with lag time
#' @examples
#' Ct <- calc_sd_1cmt_linear_oral_0_lag(t=0:24, CL=6, V=25, dur=1.5, dose=600, tlag=1.5)
#' @export
calc_sd_1cmt_linear_oral_0_lag <- function(t, dose, dur, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.24 p. 13
  Ct <- (dose / dur) * (1 / (param$k10 * param$V1)) * (1 - exp(-param$k10 * dur)) * exp(-param$k10 * (t - param$tlag - dur))
  Ct[t < param$tlag] <- 0
  Ct[t >= param$tlag & t < (dur+param$tlag)] <- (dose / dur) * (1 / (param$k10 * param$V1)) * (1 - exp(-param$k10 * (t[t >= param$tlag & t < (dur+param$tlag)] - param$tlag)))
  Ct
}

# Steady-State #####

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
#' @author Bill Denney, \email{wdenney@@humanpredictions.com}
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
calc_ss_1cmt_linear_bolus <- function(tad, tau, dose, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.1 p. 5
  Ct <- (dose/param$V1) * (exp(-param$k10*tad))/(1 - exp(-param$k10 * tau))
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with infusion at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_infusion(tad=0:36, CL=2, V=25, dose=600, tinf=1, tau=24)
#' @export
calc_ss_1cmt_linear_infusion <- function(tad, tau, dose, tinf, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.8 p. 8
  Ct <-
    (dose/tinf) *
    (1 / (param$k10*param$V1) ) *
    ((1 - exp(-param$k10 * tinf)) *
       exp(-param$k10 * (tad - tinf))) /
    (1 - exp(-param$k10 * tau))
  Ct[tad <= tinf] <-
    (dose/tinf) *
    (1 / (param$k10*param$V1) ) *
    ((1 - exp(-param$k10 * tad[tad<=tinf])) +
       exp(-param$k10 * tau) *
       (((1 - exp(-param$k10 * tinf)) *
           exp(-param$k10 * (tad[tad<=tinf] - tinf))) /
          (1 - exp(-param$k10 * tau))
       )
    )
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_0(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24)
#' @export
calc_ss_1cmt_linear_oral_0 <- function(tad, tau, dose, dur, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) at steady state - eq 1.23 p. 14
  Ct <-
    (dose / dur) *
    (1 / (param$k10 * param$V1)) *
    ((1- exp(-param$k10 * dur)) *
       exp(-param$k10 * (tad - dur)) /
       (1 - exp(-param$k10 * tau))
    )
  Ct[tad <= dur] <-
    (dose / dur) *
    (1 / (param$k10 * param$V1)) *
    ((1 - exp(-param$k10 * tad[tad <= dur])) +
       exp(-param$k10 * tau)*
       ((1 - exp(-param$k10 * dur)) *
          exp(-param$k10 * (tad[tad <= dur] - dur)) /
          (1 - exp(-param$k10 * tau)))
    )
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with zero-order oral absorption at steady state, with lag time
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_0_lag(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24, tlag=1.5)
#' @export
calc_ss_1cmt_linear_oral_0_lag <- function(tad, tau, dose, dur, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) at steady state - eq 1.26 p. 14
  Ct <-
    (dose / dur) *
    (1 / (param$k10 * param$V1)) *
    ((1 - exp(-param$k10 * dur)) *
       exp(-param$k10 * (tad - param$tlag - dur)) /
       (1 - exp(-param$k10 * tau))
    )
  Ct[tad < param$tlag] <-
    (dose / dur) *
    (1 / (param$k10 * param$V1)) *
    ((1 - exp(-param$k10 * dur)) *
       exp(-param$k10 * (tad[tad < param$tlag] + tau - param$tlag - dur)) /
       (1 - exp(-param$k10 * tau))
    )
  Ct[tad >= param$tlag & tad < (dur+param$tlag)] <-
    (dose / dur) *
    (1 / (param$k10 * param$V1)) *
    ((1 - exp(-param$k10 * (tad[tad >= param$tlag & tad < (dur+param$tlag)] - param$tlag))) +
       exp(-param$k10 * tau) *
       ((1 - exp(-param$k10 * dur)) *
          exp(-param$k10 * (tad[tad >= param$tlag & tad < (dur+param$tlag)] - param$tlag - dur)) /
          (1 - exp(-param$k10 * tau)))
    )
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state, with lag time
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_1_lag(tad=0:36, CL=2, V=25, dose=600,
#'     ka=0.25, tlag=0.75, tau=24)
#' @export
calc_ss_1cmt_linear_oral_1_lag <- function(tad, tau, dose, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.16 p. 11
  Ct <-
    (dose /param$V1) *
    (param$ka / (param$ka - param$k10)) *
    ((exp(-param$k10 * (tad - param$tlag)) / (1 - exp(-param$k10 * tau))) -
       (exp(-param$ka * (tad - param$tlag)) / (1 - exp(-param$ka * tau)))
    )
  Ct[tad < param$tlag] <-
    (dose /param$V1) *
    (param$ka / (param$ka - param$k10)) *
    ((exp(-param$k10 * (tad[tad < param$tlag] + tau - param$tlag)) /
        (1 - exp(-param$k10 * tau))) -
       (exp(-param$ka * (tad[tad < param$tlag] + tau - param$tlag)) /
          (1 - exp(-param$ka * tau)))
    )
  Ct
}

#' @describeIn calc_ss_1cmt
#' Calculate C(t) for a 1-compartment linear model with first-order oral absorption at steady state
#' @examples
#' Ct <- calc_ss_1cmt_linear_oral_1(tad=0:36, CL=2, V=25, dose=600, ka=0.25, tau=24)
#' @export
calc_ss_1cmt_linear_oral_1 <- function(tad, tau, dose, ...) {
  param <- calc_derived_1cpt(..., sigdig=Inf)
  ### C(t) after single dose - eq 1.8 p. 8
  Ct <-
    (dose/param$V1) *
    (param$ka / (param$ka - param$k10)) *
    ((exp(-param$k10 * tad) /
        (1 - exp(-param$k10 * tau))) -
       (exp(-param$ka * tad) / (1 - exp(-param$ka * tau)))
    )
  Ct
}
