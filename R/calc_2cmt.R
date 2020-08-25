# Single-Dose ####

#' Calculate C(t) for a 1-compartment linear model
#'
#' @inheritParams calc_sd_1cmt
#' @param ... Passed to `calc_derived_2cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_sd_2cmt <- function(t, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single IV bolus dose
#' @examples
#' Ct <- calc_sd_2cmt_linear_bolus(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10)
#' @export
calc_sd_2cmt_linear_bolus <- function(t, dose, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.2 p. 21
  A <- (1/param$V1) * ((param$alpha - param$k21)/(param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21)/(param$beta - param$alpha))
  ### C(t) after single dose - eq 1.36 p. 22
  Ct <- dose * (A * exp(-param$alpha * t) + (B * exp(-param$beta * t)))
  Ct
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single first-order oral dose with a lag time
#' @examples
#' Ct <- calc_sd_2cmt_linear_oral_1_lag(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tlag = 2)
#' @export
calc_sd_2cmt_linear_oral_1_lag <- function(t, dose, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.3 p. 24
  A <- (param$ka/param$V1) * ((param$k21 - param$alpha) / ((param$ka - param$alpha) * (param$beta - param$alpha)))
  B <- (param$ka/param$V1) * ((param$k21 - param$beta) / ((param$ka - param$beta) * (param$alpha - param$beta)))
  ### C(t) after single dose - eq 1.41 p. 25
  Ct <- dose * ((A * exp(-param$alpha * (t - param$tlag))) + (B * exp(-param$beta * (t - param$tlag))) - ((A + B) * exp(-param$ka * (t - param$tlag))))
  Ct[t < param$tlag] <- 0
  Ct
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single infusion
#' @examples
#' Ctrough <- calc_sd_2cmt_linear_infusion(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10, tinf = 1)
#' @export
calc_sd_2cmt_linear_infusion <- function(t, dose, tinf, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.2 p. 21
  A <- (1/param$V1) * ((param$alpha - param$k21)/(param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21)/(param$beta - param$alpha))
  ### C(t) after single dose - eq 1.36 p. 22
  Ct <- (dose/tinf) * (((A/param$alpha) * (((1 - exp(-param$alpha * tinf)) * (exp(-param$alpha * (t-tinf)))))) +
                         ((B/param$beta) * (((1 - exp(-param$beta * tinf)) * (exp(-param$beta * (t-tinf)))))))
  Ct[t <= tinf] <- (dose/tinf) * (((A/param$alpha) * (1 - exp(-param$alpha * t[t <= tinf]))) +
                                    ((B/param$beta) * (1 - exp(-param$beta * t[t <= tinf]))))
  Ct
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single zero-order oral dose, with lag time
#' @examples
#' Ctrough <- calc_sd_2cmt_linear_oral_0_lag(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tlag=2)
#' @export
calc_sd_2cmt_linear_oral_0_lag <- function(t, dose, dur, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.4 p. 28
  A <- (1/param$V1) * ((param$alpha - param$k21) / (param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21) / (param$beta - param$alpha))
  ### C(t) after single dose - eq 1.54 p. 31
  Ct <- (dose / dur) * (((A / param$alpha) * (1 - exp(-param$alpha * dur)) * exp(-param$alpha * (t - param$tlag - dur))) +
                          ((B / param$beta) * (1 - exp(-param$beta * dur)) * exp(-param$beta * (t - param$tlag - dur))))
  Ct[t < param$tlag] <- 0
  Ct[t >= param$tlag & t < (dur+param$tlag)] <- (dose / dur) * ((A / param$alpha) * (1 - exp(-param$alpha * (t[t >= param$tlag & t < (dur+param$tlag)] - param$tlag))) +
                                                     (B / param$beta) * (1 - exp(-param$beta * (t[t >= param$tlag & t < (dur+param$tlag)] - param$tlag))))
  Ct
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single zero-order oral dose, with lag time
#' @examples
#' Ct <- calc_sd_2cmt_linear_oral_0_lag(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tlag=2)
#' @export
calc_sd_2cmt_linear_oral_0_lag <- function(t, dose, dur, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.4 p. 28
  A <- (1/param$V1) * ((param$alpha - param$k21) / (param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21) / (param$beta - param$alpha))
  ### C(t) after single dose - eq 1.54 p. 31
  Ct <- (dose / dur) * (((A / param$alpha) * (1 - exp(-param$alpha * dur)) * exp(-param$alpha * (t - param$tlag - dur))) +
                          ((B / param$beta) * (1 - exp(-param$beta * dur)) * exp(-param$beta * (t - param$tlag - dur))))
  Ct[t < param$tlag] <- 0
  Ct[t >= param$tlag & t < (dur+param$tlag)] <- (dose / dur) * ((A / param$alpha) * (1 - exp(-param$alpha * (t[t >= param$tlag & t < (dur+param$tlag)] - param$tlag))) +
                                                     (B / param$beta) * (1 - exp(-param$beta * (t[t >= param$tlag & t < (dur+param$tlag)] - param$tlag))))
  Ct
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single first-order oral dose
#' @examples
#' Ct <- calc_sd_2cmt_linear_oral_1(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1)
#' @export
calc_sd_2cmt_linear_oral_1 <- function(t, dose, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.3 p. 24
  A <- (param$ka/param$V1) * ((param$k21 - param$alpha) / ((param$ka - param$alpha) * (param$beta - param$alpha)))
  B <- (param$ka/param$V1) * ((param$k21 - param$beta) / ((param$ka - param$beta) * (param$alpha - param$beta)))
  ### C(t) after single dose - eq 1.41 p. 25
  Ct <- dose * ((A * exp(-param$alpha * t)) + (B * exp(-param$beta * t)) - ((A + B) * exp(-param$ka * t)))
  Ct
}

#' @describeIn calc_sd_2cmt
#' Calculate C(t) for a 2-compartment linear model after a single zero-order oral dose
#' @examples
#' Ct <- calc_sd_2cmt_linear_oral_0(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1)
#' @export
calc_sd_2cmt_linear_oral_0 <- function(t, dose, dur, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.4 p. 28
  A <- (1/param$V1) * ((param$alpha - param$k21) / (param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21) / (param$beta - param$alpha))
  ### C(t) after single dose - eq 1.51 p. 29
  Ct <- (dose / dur) * (((A / param$alpha) * (1 - exp(-param$alpha * dur)) * exp(-param$alpha * (t - dur))) +
                          ((B / param$beta) * (1 - exp(-param$beta * dur)) * exp(-param$beta * (t - dur))))
  Ct[t < dur] <- (dose / dur) * ((A / param$alpha) * (1 - exp(-param$alpha * t[t < dur])) +
                                   (B / param$beta) * (1 - exp(-param$beta * t[t < dur])))
  Ct
}

# Steady-state ####
#' Calculate C(t) for a 2-compartment linear model at steady-state
#'
#' @inheritParams calc_ss_1cmt
#' @param ... Passed to `calc_derived_2cpt()`
#' @return Concentration of drug at requested time (\code{t}) at steady-state, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_ss_2cmt <- function(tad, tau, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model with IV bolus dosing at steady-state
#' @examples
#' Ct <- calc_ss_2cmt_linear_bolus(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10, tau=24)
#' @export
calc_ss_2cmt_linear_bolus <- function(tad, tau, dose, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.2 p. 21
  A <- (1/param$V1) * ((param$alpha - param$k21)/(param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21)/(param$beta - param$alpha))
  ### C(t) after single dose - eq 1.36 p. 22
  Ct <-
    dose *
    ((A * exp(-param$alpha * tad) /
        (1 - exp(-param$alpha * tau))) +
       ((B * exp(-param$beta * tad) /
           (1 - exp(-param$beta * tau))))
    )
  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model with infusion at steady state
#' @examples
#' Ct <- calc_ss_2cmt_linear_infusion(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10, tinf = 1, tau = 12)
#' @export
calc_ss_2cmt_linear_infusion <- function(tad, tau, dose, tinf, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.2 p. 21
  A <- 1 / param$V1 * ((param$alpha-param$k21)/(param$alpha-param$beta))
  B <- 1 / param$V1 * ((param$beta-param$k21)/(param$beta-param$alpha))
  ### C(t) at steady state - eq 1.38 p. 23
  InfA1 <- 1-exp(-1*param$alpha*tad)
  InfA2 <- exp(-1*param$alpha*tau)*(((1-exp(-1*param$alpha*tinf))* exp(-1*param$alpha*(tad-tinf)))/(1-exp(-1*param$alpha*tau)))
  InfB1 <- 1-exp(-1*param$beta*tad)
  InfB2 <- exp(-1*param$beta*tau)*(((1-exp(-1*param$beta*tinf))* exp(-1*param$beta*(tad-tinf)))/(1-exp(-1*param$beta*tau)))
  
  PInfA2 <- (((1-exp(-1*param$alpha*tinf))* exp(-1*param$alpha*(tad-tinf)))/(1-exp(-1*param$alpha*tau)))
  PInfB2 <- (((1-exp(-1*param$beta *tinf))* exp(-1*param$beta *(tad-tinf)))/(1-exp(-1*param$beta *tau)))
  
  InfProf  <- A/param$alpha*(InfA1+InfA2)+B/param$beta*(InfB1+InfB2)
  PInfProf <- A/param$alpha*(PInfA2)+B/param$beta*(PInfB2)
  
  Ct <- dose/tinf*(InfProf*(tad<=tinf)+PInfProf*(tad>tinf))
  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with zero-order oral dosing
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_0(tad = 23, CL = 2.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tau = 24)
#' @export
calc_ss_2cmt_linear_oral_0 <- function(tad, tau, dose, dur, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.4 p. 28
  A <- (1/param$V1) * ((param$alpha - param$k21) / (param$alpha - param$beta ))
  B <- (1/param$V1) * ((param$beta  - param$k21) / (param$beta  - param$alpha))
  ### C(t) at steady state - eq 1.53 p. 30
  c1  <- dose / dur
  ca1 <- A / param$alpha
  cb1 <- B / param$beta
  
  c1a1 <- 1 - exp(-param$alpha * dur)
  c1a2 <- exp(-param$alpha * (tad - dur))
  c1a3 <- 1 - exp(-param$alpha * tau)
  
  c1b1 <- 1 - exp(-param$beta * dur)
  c1b2 <- exp(-param$beta * (tad - dur))
  c1b3 <- 1 - exp(-param$beta * tau)
  
  c2a1 <- 1 - exp(-param$alpha * tad[tad < dur])
  c2a2 <- exp(-param$alpha * tau)
  c2a3 <- 1 - exp(-param$alpha * dur)
  c2a4 <- exp(-param$alpha * (tad[tad < dur] - dur))
  c2a5 <- 1 - exp(-param$alpha * tau)
  
  c2b1 <- 1 - exp(-param$beta * tad[tad < dur])
  c2b2 <- exp(-param$beta * tau)
  c2b3 <- 1 - exp(-param$beta * dur)
  c2b4 <- exp(-param$beta * (tad[tad < dur] - dur))
  c2b5 <- 1 - exp(-param$beta * tau)
  
  Ct <- c1 * ((ca1 * ((c1a1 * c1a2) / c1a3)) +
                (cb1 * ((c1b1 * c1b2) / c1b3)))
  Ct[tad < dur] <- c1 * ((ca1 * (c2a1 + c2a2*((c2a3 * c2a4)/c2a5))) +
                           (cb1 * (c2b1 + c2b2*((c2b3 * c2b4)/c2b5))))
  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with first-order oral dosing
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_1_lag(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tau=24, tlag=2)
#' @export
calc_ss_2cmt_linear_oral_1_lag <- function(tad, tau, dose, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.3 p. 24
  A <- (param$ka/param$V1) * ((param$k21 - param$alpha) / ((param$ka - param$alpha) * (param$beta - param$alpha)))
  B <- (param$ka/param$V1) * ((param$k21 - param$beta) / ((param$ka - param$beta) * (param$alpha - param$beta)))
  ### C(t) after single dose - eq 1.46 p. 26
  Ct <- dose * (((A * exp(-param$alpha * (tad - param$tlag))) / (1 - exp(-param$alpha * tau)))  +
                  ((B * exp(-param$beta * (tad - param$tlag))) / (1 - exp(-param$beta * tau))) -
                  (((A + B) * exp(-param$ka * (tad - param$tlag))) / (1 - exp(-param$ka * tau))))
  Ct[tad < param$tlag] <- dose * (((A * exp(-param$alpha * (tad[tad < param$tlag] + tau - param$tlag))) / (1 - exp(-param$alpha * tau)))  +
                              ((B * exp(-param$beta * (tad[tad < param$tlag] + tau - param$tlag))) / (1 - exp(-param$beta * tau))) -
                              (((A + B) * exp(-param$ka * (tad[tad < param$tlag] + tau - param$tlag))) / (1 - exp(-param$ka * tau))))
  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with zero-order oral dosing and a lag time
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_0_lag(tad = 23, CL = 2.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tau = 24, tlag=2)
#' @export
calc_ss_2cmt_linear_oral_0_lag <- function(tad, tau, dose, dur, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.4 p. 28
  A <- (1/param$V1) * ((param$alpha - param$k21) / (param$alpha - param$beta))
  B <- (1/param$V1) * ((param$beta - param$k21) / (param$beta - param$alpha))
  ### C(t) at steady state - eq 1.56 p. 33
  c1  <- dose / dur
  ca1 <- A / param$alpha
  cb1 <- B / param$beta
  
  # rest
  
  c1a1 <- 1 - exp(-param$alpha * dur)
  c1a2 <- exp(-param$alpha * (tad - param$tlag - dur))
  c1a3 <- 1 - exp(-param$alpha * tau)
  
  c1b1 <- 1 - exp(-param$beta * dur) 
  c1b2 <- exp(-param$beta * (tad - param$tlag - dur))
  c1b3 <- 1 - exp(-param$beta * tau)
  
  Ct <- c1 * (ca1 * ((c1a1 * c1a2) / c1a3) +
                cb1 * ((c1b1 * c1b2) / c1b3))
  
  # during lag
  
  c2a1 <- 1 - exp(-param$alpha * dur)
  c2a2 <- exp(-param$alpha * (tad[tad <= param$tlag] + tau - param$tlag - dur))
  c2a3 <- 1 - exp(-param$alpha * tau)
  
  c2b1 <- 1 - exp(-param$beta * dur)
  c2b2 <- exp(-param$beta * (tad[tad <= param$tlag] + tau - param$tlag - dur))
  c2b3 <- 1 - exp(-param$beta * tau)
  
  Ct[tad <= param$tlag] <- c1 * (ca1 * ((c2a1 * c2a2) / c2a3) +
                             cb1 * ((c2b1 * c2b2) / c2b3))
  
  # during zero order absorption
  
  c3a1 <- 1 - exp(-param$alpha * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag))
  c3a2 <- 1 - exp(-param$alpha * dur)
  c3a3 <- exp(-param$alpha * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag - dur))
  c3a4 <- 1 - exp(-param$alpha * tau)
  
  c3b1 <- 1 - exp(-param$beta * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag))
  c3b2 <- 1 - exp(-param$beta * dur)
  c3b3 <- exp(-param$beta * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag - dur))
  c3b4 <- 1 - exp(-param$beta * tau)
  
  Ct[tad > param$tlag & tad <= param$tlag + dur] <-
    c1 * (ca1 * (c3a1 + exp(-param$alpha * tau) * (c3a2 * c3a3 / c3a4)) +
            cb1 * (c3b1 + exp(-param$beta * tau) * (c3b2 * c3b3 / c3b4)))
  
  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with first-order oral dosing
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_1(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tau=24)
#' @export
calc_ss_2cmt_linear_oral_1 <- function(tad, tau, dose, ...) {
  param <- calc_derived_2cpt(..., sigdig=Inf)
  ### macroconstants - 1.2.3 p. 24
  A <- (param$ka/param$V1) * ((param$k21 - param$alpha) / ((param$ka - param$alpha) * (param$beta - param$alpha)))
  B <- (param$ka/param$V1) * ((param$k21 - param$beta) / ((param$ka - param$beta) * (param$alpha - param$beta)))
  ### C(t) after single dose - eq 1.43 p. 25
  Ct <- dose * (((A * exp(-param$alpha * tad)) / (1 - exp(-param$alpha * tau)))  + ((B * exp(-param$beta * tad)) / (1 - exp(-param$beta * tau))) - (((A + B) * exp(-param$ka * tad)) / (1 - exp(-param$ka * tau))))
  Ct
}
