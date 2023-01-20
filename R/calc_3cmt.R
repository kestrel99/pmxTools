# Single Dose ####

#' Calculate C(t) for a 3-compartment linear model
#'
#' @inheritParams calc_sd_1cmt
#' @param ... Passed to `calc_derived_3cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @author Bill Denney, \email{wdenney@@humanpredictions.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_sd_3cmt <- function(t, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_sd_3cmt
#' Calculate C(t) for a 3-compartment linear model after a single IV bolus dose
#' @examples
#' Ct <- calc_sd_3cmt_linear_bolus(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dose = 100)
#' @export
calc_sd_3cmt_linear_bolus <- function(t, dose, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.1 p. 39
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.61 p. 40
  Ct <- dose * ((A * exp(-param$alpha * t)) + (B * exp(-param$beta * t)) + (C * exp(-param$gamma * t)))
  Ct
}

#' @describeIn calc_sd_3cmt
#' Calculate C(t) for a 3-compartment linear model after a single oral dose
#' @examples
#' Ct <- calc_sd_3cmt_linear_oral_1_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tlag = 1.5)
#' @export
calc_sd_3cmt_linear_oral_1_lag <- function(t, dose, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.3 p. 44
  A <- (1/param$V1) * (param$ka / (param$ka - param$alpha)) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * (param$ka / (param$ka - param$beta))  * ((param$k21 - param$beta)/(param$beta - param$alpha))  * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * (param$ka / (param$ka - param$gamma)) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose with lag time - eq 1.74 p. 45
  Ct <- dose * ((A * exp(-param$alpha*(t-param$tlag))) + (B * exp(-param$beta*(t-param$tlag))) + (C * exp(-param$gamma*(t-param$tlag))) - ((A + B + C)*exp(-param$ka*(t-param$tlag))))
  Ct[t<param$tlag] <- 0
  Ct
}

#' @describeIn calc_sd_3cmt
#' Calculate C(t) for a 3-compartment linear model after a single IV infusion
#' @examples
#' Ct <- calc_sd_3cmt_linear_infusion(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dose = 100, tinf=1)
#' @export
calc_sd_3cmt_linear_infusion <- function(t, dose, tinf, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.2 p. 41
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.66 p. 42
  
  c1 <- dose / tinf
  ca <- A / param$alpha
  cb <- B / param$beta
  cc <- C / param$gamma
  
  c1a1 <- 1 - exp(-param$alpha * tinf)
  c1a2 <- exp(-param$alpha * (t - tinf))
  c1b1 <- 1 - exp(-param$beta * tinf)
  c1b2 <- exp(-param$beta * (t - tinf))
  c1c1 <- 1 - exp(-param$gamma * tinf)
  c1c2 <- exp(-param$gamma * (t - tinf))
  
  Ct <- c1 * ((ca * c1a1 * c1a2) + (cb * c1b1 * c1b2) + (cc * c1c1 * c1c2))
  
  c2a1 <- 1 - exp(-param$alpha * t[t <= tinf])
  c2b1 <- 1 - exp(-param$beta * t[t <= tinf])
  c2c1 <- 1 - exp(-param$gamma * t[t <= tinf])
  
  Ct[t <= tinf] <- c1 * ((ca * c2a1) + (cb * c2b1) + (cc * c2c1))
  
  Ct
}

#' @describeIn calc_sd_3cmt
#' Calculate C(t) for a 3-compartment linear model after a single dose, with zero-order absorption
#' @examples
#' Ct <- calc_sd_3cmt_linear_oral_0(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100)
#' @export
calc_sd_3cmt_linear_oral_0 <- function(t, dose, dur, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.4 p. 47
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.81 p. 48
  Ct <-
    (dose / dur) *
    ((A/param$alpha) * (1 - exp(-param$alpha * dur)) * exp(-param$alpha * (t - dur)) +
       (B/param$beta) * (1 - exp(-param$beta * dur)) * exp(-param$beta * (t - dur)) +
       (C/param$gamma) * (1 - exp(-param$gamma * dur)) * exp(-param$gamma * (t - dur)))
  Ct[t < dur] <-
    (dose / dur) *
    ((A/param$alpha) * (1 - exp(-param$alpha * t[t < dur])) +
       (B/param$beta) * (1 - exp(-param$beta * t[t < dur])) +
       (C/param$gamma) * (1 - exp(-param$gamma * t[t < dur])))
  Ct
}

#' @describeIn calc_sd_3cmt
#' Calculate C(t) for a 3-compartment linear model after a single dose, with zero-order absorption and a lag time
#' @examples
#' Ct <- calc_sd_3cmt_linear_oral_0_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tlag=1.5)
#' @export
calc_sd_3cmt_linear_oral_0_lag <- function(t, dose, dur, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.4 p. 47
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.84 p. 49
  Ct <-
    (dose / dur) *
    ((A/param$alpha) * (1 - exp(-param$alpha * dur)) * exp(-param$alpha * (t - param$tlag - dur)) +
       (B/param$beta) *  (1 - exp(-param$beta * dur))  * exp(-param$beta  * (t - param$tlag - dur)) +
       (C/param$gamma) * (1 - exp(-param$gamma * dur)) * exp(-param$gamma * (t - param$tlag - dur)))
  Ct[t <= dur + param$tlag & t > param$tlag] <-
    (dose / dur) *
    ((A/param$alpha) * (1 - exp(-param$alpha * (t[t <= dur + param$tlag & t > param$tlag] - param$tlag))) +
       (B/param$beta)  * (1 - exp(-param$beta  * (t[t <= dur + param$tlag & t > param$tlag] - param$tlag))) +
       (C/param$gamma) * (1 - exp(-param$gamma * (t[t <= dur + param$tlag & t > param$tlag] - param$tlag))))
  Ct[t <= param$tlag] <- 0
  Ct
}

#' @describeIn calc_sd_3cmt
#' Calculate C(t) for a 3-compartment linear model after a single oral dose
#' @examples
#' Ct <- calc_sd_3cmt_linear_oral_1(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100)
#' @export
calc_sd_3cmt_linear_oral_1 <- function(t, dose, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.3 p. 44
  A <- (1/param$V1) * (param$ka / (param$ka - param$alpha)) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * (param$ka / (param$ka - param$beta))  * ((param$k21 - param$beta)/(param$beta - param$alpha))  * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * (param$ka / (param$ka - param$gamma)) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.71 p. 44
  Ct <-
    dose *
    ((A * exp(-param$alpha*t)) +
       (B * exp(-param$beta*t)) +
       (C * exp(-param$gamma*t)) -
       ((A + B + C)*exp(-param$ka*t))
    )
  Ct
}

# Steady-State ####

#' Calculate C(t) for a 3-compartment linear model at steady-state
#'
#' @inheritParams calc_ss_1cmt
#' @param ... Passed to `calc_derived_3cpt()`
#' @return Concentration of drug at requested time (\code{t}) at steady-state, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_ss_3cmt <- function(tad, tau, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady state with IV bolus dosing
#' @examples
#' Ct <- calc_ss_3cmt_linear_bolus(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dose = 100, tau=24)
#' @export
calc_ss_3cmt_linear_bolus <- function(tad, tau, dose, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.1 p. 39
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.61 p. 40
  Ct <-
    dose *
    ((A * (exp(-param$alpha * tad) / (1 - exp(-param$alpha * tau)))) +
       (B * (exp(-param$beta * tad) / (1 - exp(-param$beta * tau)))) +
       (C * (exp(-param$gamma * tad) / (1 - exp(-param$gamma * tau)))))
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady-state with first-order oral dosing with a lag time
#' @examples
#' Ctrough <- calc_ss_3cmt_linear_oral_1_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau=24, tlag = 1.5)
#' @export
calc_ss_3cmt_linear_oral_1_lag <- function(tad, tau, dose, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.3 p. 44
  A <- (1/param$V1) * (param$ka / (param$ka - param$alpha)) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * (param$ka / (param$ka - param$beta))  * ((param$k21 - param$beta)/(param$beta - param$alpha))  * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * (param$ka / (param$ka - param$gamma)) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose with lag time - eq 1.76 p. 45
  Ct <-
    dose *
    (((A * exp(-param$alpha*(tad-param$tlag)))/(1-exp(-param$alpha*tau))) +
       ((B * exp(-param$beta*(tad-param$tlag)))/(1-exp(-param$beta*tau))) +
       ((C * exp(-param$gamma*(tad-param$tlag)))/(1-exp(-param$gamma*tau))) -
       (((A + B + C)*exp(-param$ka*(tad-param$tlag)))/(1-exp(-param$ka*tau))))
  Ct[tad<param$tlag] <-
    dose *
    (((A * exp(-param$alpha*(tad[tad<param$tlag]+tau-param$tlag)))/(1-exp(-param$alpha*tau))) +
       ((B * exp(-param$beta*(tad[tad<param$tlag]+tau-param$tlag)))/(1-exp(-param$beta*tau))) +
       ((C * exp(-param$gamma*(tad[tad<param$tlag]+tau-param$tlag)))/(1-exp(-param$gamma*tau))) -
       (((A + B + C)*exp(-param$ka*(tad[tad<param$tlag]+tau-param$tlag)))/(1-exp(-param$ka*tau))))
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady state with IV infusions
#' @examples
#' Ct <- calc_ss_3cmt_linear_infusion(tad = 11.75, CL = 2.5, V1 = 20, V2 = 50,
#'     V3 = 100, Q2 = 0.5, Q3 = 0.05, dose = 1000, tinf=1, tau=24)
#' @export
calc_ss_3cmt_linear_infusion <- function(tad, tau, dose, tinf, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.2 p. 41
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.68 p. 43
  c1 <- dose / tinf
  ca <- A / param$alpha
  cb <- B / param$beta
  cc <- C / param$gamma
  
  c1a1 <- 1 - exp(-param$alpha * tinf)
  c1a2 <- exp(-param$alpha * (tad - tinf))
  c1a3 <- 1 - exp(-param$alpha * tau)
  
  c1b1 <- 1 - exp(-param$beta * tinf)
  c1b2 <- exp(-param$beta * (tad - tinf))
  c1b3 <- 1 - exp(-param$beta * tau)
  
  c1c1 <- 1 - exp(-param$gamma * tinf)
  c1c2 <- exp(-param$gamma * (tad - tinf))
  c1c3 <- 1 - exp(-param$gamma * tau)
  
  Ct <-
    c1 *
    ((ca * ((c1a1 * c1a2)/c1a3)) +
       (cb * ((c1b1 * c1b2)/c1b3)) +
       (cc * ((c1c1 * c1c2)/c1c3)))
  
  c2a1 <- 1 - exp(-param$alpha * tad[tad <= tinf])
  c2a2 <- exp(-param$alpha * tau)
  c2a3 <- 1 - exp(-param$alpha * tinf)
  c2a4 <- exp(-param$alpha * (tad[tad <= tinf] - tinf))
  c2a5 <- 1 - exp(-param$alpha * tau)
  
  c2b1 <- 1 - exp(-param$beta * tad[tad <= tinf])
  c2b2 <- exp(-param$beta * tau)
  c2b3 <- 1 - exp(-param$beta * tinf)
  c2b4 <- exp(-param$beta * (tad[tad <= tinf] - tinf))
  c2b5 <- 1 - exp(-param$beta * tau)
  
  c2c1 <- 1 - exp(-param$gamma * tad[tad <= tinf])
  c2c2 <- exp(-param$gamma * tau)
  c2c3 <- 1 - exp(-param$gamma * tinf)
  c2c4 <- exp(-param$gamma * (tad[tad <= tinf] - tinf))
  c2c5 <- 1 - exp(-param$gamma * tau)
  
  Ct[tad <= tinf] <-
    c1 *
    (ca * (c2a1 + (c2a2 * ((c2a3 * c2a4) / c2a5))) +
       cb * (c2b1 + (c2b2 * ((c2b3 * c2b4) / c2b5))) +
       cc * (c2c1 + (c2c2 * ((c2c3 * c2c4) / c2c5))))
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady state, with zero-order absorption
#' @examples
#' Ct <- calc_ss_3cmt_linear_oral_0(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24)
#' @export
calc_ss_3cmt_linear_oral_0 <- function(tad, tau, dose, dur, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.4 p. 47
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) after single dose - eq 1.83 p. 49
  Ct <-
    (dose / dur) *
    ((A/param$alpha)   * ((1 - exp(-param$alpha * dur)) * exp(-param$alpha * (tad - dur)))/(1 - exp(-param$alpha * tau)) +
       (B/param$beta)  * ((1 - exp(-param$beta  * dur)) * exp(-param$beta  * (tad - dur)))/(1 - exp(-param$beta  * tau)) +
       (C/param$gamma) * ((1 - exp(-param$gamma * dur)) * exp(-param$gamma * (tad - dur)))/(1 - exp(-param$gamma * tau)))
  
  tm1a <- 1 - exp(-param$alpha * tad[tad <= dur])
  tm2a <- exp(-param$alpha * tau)
  tm3a <- 1 - exp(-param$alpha * dur)
  tm4a <- exp(-param$alpha * (tad[tad <= dur] - dur))
  tm5a <- 1 - exp(-param$alpha * tau)
  
  tm1b <- 1 - exp(-param$beta * tad[tad <= dur])
  tm2b <- exp(-param$beta * tau)
  tm3b <- 1 - exp(-param$beta * dur)
  tm4b <- exp(-param$beta * (tad[tad <= dur] - dur))
  tm5b <- 1 - exp(-param$beta * tau)
  
  tm1c <- 1 - exp(-param$gamma * tad[tad <= dur])
  tm2c <- exp(-param$gamma * tau)
  tm3c <- 1 - exp(-param$gamma * dur)
  tm4c <- exp(-param$gamma * (tad[tad <= dur] - dur))
  tm5c <- 1 - exp(-param$gamma * tau)
  
  Ct[tad <= dur] <-
    (dose / dur) *
    ((A/param$alpha) * (tm1a + (tm2a * ((tm3a * tm4a)/tm5a))) +
       (B/param$beta)  * (tm1b + (tm2b * ((tm3b * tm4b)/tm5b))) +
       (C/param$gamma) * (tm1c + (tm2c * ((tm3c * tm4c)/tm5c))))
  
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady state, with zero-order absorption and lag time
#' @examples
#' Ct <- calc_ss_3cmt_linear_oral_0_lag(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24, tlag = 1.5)
#' @export
calc_ss_3cmt_linear_oral_0_lag <- function(tad, tau, dose, dur, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.4 p. 47
  A <- (1/param$V1) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * ((param$k21 - param$beta)/(param$beta - param$alpha)) * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) at steady state - eq 1.86 p. 51
  
  tm1a1 <- 1 - exp(-param$alpha * dur)
  tm2a1 <- exp(-param$alpha * (tad - param$tlag - dur))
  tm3a1 <- 1 - exp(-param$alpha * tau)
  
  tm1b1 <- 1 - exp(-param$beta * dur)
  tm2b1 <- exp(-param$beta * (tad - param$tlag - dur))
  tm3b1 <- 1 - exp(-param$beta * tau)
  
  tm1c1 <- 1 - exp(-param$gamma * dur)
  tm2c1 <- exp(-param$gamma * (tad - param$tlag - dur))
  tm3c1 <- 1 - exp(-param$gamma * tau)
  
  Ct <-
    (dose / dur) *
    ((A/param$alpha) * ((tm1a1 * tm2a1)/tm3a1) +
       (B/param$beta)  * ((tm1b1 * tm2b1)/tm3b1) +
       (C/param$gamma) * ((tm1c1 * tm2c1)/tm3c1))
  
  
  tm1a2 <- 1 - exp(-param$alpha * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag))
  tm2a2 <- exp(-param$alpha * tau)
  tm3a2 <- 1 - exp(-param$alpha * dur)
  tm4a2 <- exp(-param$alpha * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag - dur))
  tm5a2 <- 1 - exp(-param$alpha * tau)
  
  tm1b2 <- 1 - exp(-param$beta * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag))
  tm2b2 <- exp(-param$beta * tau)
  tm3b2 <- 1 - exp(-param$beta * dur)
  tm4b2 <- exp(-param$beta * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag - dur))
  tm5b2 <- 1 - exp(-param$beta * tau)
  
  tm1c2 <- 1 - exp(-param$gamma * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag))
  tm2c2 <- exp(-param$gamma * tau)
  tm3c2 <- 1 - exp(-param$gamma * dur)
  tm4c2 <- exp(-param$gamma * (tad[tad > param$tlag & tad <= param$tlag + dur] - param$tlag - dur))
  tm5c2 <- 1 - exp(-param$gamma * tau)
  
  Ct[tad > param$tlag & tad <= param$tlag + dur] <-
    (dose/dur) *
    ((A/param$alpha) * (tm1a2 + (tm2a2 * (tm3a2 * tm4a2)/tm5a2)) +
       (B/param$beta)  * (tm1b2 + (tm2b2 * (tm3b2 * tm4b2)/tm5b2)) +
       (C/param$gamma) * (tm1c2 + (tm2c2 * (tm3c2 * tm4c2)/tm5c2)))
  
  tm1a3 <- 1 - exp(-param$alpha * dur)
  tm2a3 <- exp(-param$alpha * (tad[tad <= param$tlag] + tau - param$tlag - dur))
  tm3a3 <- 1 - exp(-param$alpha * tau)
  
  tm1b3 <- 1 - exp(-param$beta * dur)
  tm2b3 <- exp(-param$beta * (tad[tad <= param$tlag] + tau - param$tlag - dur))
  tm3b3 <- 1 - exp(-param$beta * tau)
  
  tm1c3 <- 1 - exp(-param$gamma * dur)
  tm2c3 <- exp(-param$gamma * (tad[tad <= param$tlag] + tau - param$tlag - dur))
  tm3c3 <- 1 - exp(-param$gamma * tau)
  
  Ct[tad <= param$tlag] <-
    (dose/dur) *
    ((A/param$alpha) * ((tm1a3 * tm2a3) / tm3a3) +
       (B/param$beta)  * ((tm1b3 * tm2b3) / tm3b3) +
       (C/param$gamma) * ((tm1c3 * tm2c3) / tm3c3))
  
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady-state with first-order oral dosing
#' @examples
#' Ct <- calc_ss_3cmt_linear_oral_1(tad = 11.75, CL = 3.5, V1 = 20,
#'     V2 = 500, V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau = 24)
#' @export
calc_ss_3cmt_linear_oral_1 <- function(tad, tau, dose, ...) {
  param <- calc_derived_3cpt(..., sigdig=Inf)
  ### macroconstants - 1.3.3 p. 44
  A <- (1/param$V1) * (param$ka / (param$ka - param$alpha)) * ((param$k21 - param$alpha)/(param$alpha - param$beta)) * ((param$k31 - param$alpha)/(param$alpha - param$gamma))
  B <- (1/param$V1) * (param$ka / (param$ka - param$beta))  * ((param$k21 - param$beta)/(param$beta - param$alpha))  * ((param$k31 - param$beta)/(param$beta - param$gamma))
  C <- (1/param$V1) * (param$ka / (param$ka - param$gamma)) * ((param$k21 - param$gamma)/(param$gamma - param$beta)) * ((param$k31 - param$gamma)/(param$gamma - param$alpha))
  ### C(t) at steady state - eq 1.73 p. 45
  Ct <-
    dose *
    ((A * exp(-param$alpha*tad))/(1-exp(-param$alpha*tau)) +
       (B * exp(-param$beta*tad))/(1-exp(-param$beta*tau)) +
       (C * exp(-param$gamma*tad))/(1-exp(-param$gamma*tau)) -
       ((A + B + C)*exp(-param$ka*tad))/(1-exp(-param$ka*tau)))
  Ct
}
