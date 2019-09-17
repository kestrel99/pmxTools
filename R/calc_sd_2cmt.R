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
  Ct[t >= param$tlag & t < dur] <- (dose / dur) * ((A / param$alpha) * (1 - exp(-param$alpha * (t[t >= param$tlag & t < dur] - param$tlag))) +
                                                     (B / param$beta) * (1 - exp(-param$beta * (t[t >= param$tlag & t < dur] - param$tlag))))
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
  Ct[t >= param$tlag & t < dur] <- (dose / dur) * ((A / param$alpha) * (1 - exp(-param$alpha * (t[t >= param$tlag & t < dur] - param$tlag))) +
                                                     (B / param$beta) * (1 - exp(-param$beta * (t[t >= param$tlag & t < dur] - param$tlag))))
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
