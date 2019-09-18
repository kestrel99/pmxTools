#' Calculate C(t) for a 1-compartment linear model
#'
#' @inheritParams calc_sd_1cmt
#' @param ... Passed to `calc_derived_3cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
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
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.61 p. 40
  Ct <- dose * ((A * exp(-alpha * t)) + (B * exp(-beta * t)) + (C * exp(-gamma * t)))
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
  A <- (1/V1) * (ka / (ka - alpha)) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * (ka / (ka - beta))  * ((k21 - beta)/(beta - alpha))  * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * (ka / (ka - gamma)) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose with lag time - eq 1.74 p. 45
  Ct <- dose * ((A * exp(-alpha*(t-tlag))) + (B * exp(-beta*(t-tlag))) + (C * exp(-gamma*(t-tlag))) - ((A + B + C)*exp(-ka*(t-tlag))))
  Ct[t<tlag] <- 0
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
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.66 p. 42
  
  c1 <- dose / tinf
  ca <- A / alpha
  cb <- B / beta
  cc <- C / gamma
  
  c1a1 <- 1 - exp(-alpha * tinf)
  c1a2 <- exp(-alpha * (t - tinf))
  c1b1 <- 1 - exp(-beta * tinf)
  c1b2 <- exp(-beta * (t - tinf))
  c1c1 <- 1 - exp(-gamma * tinf)
  c1c2 <- exp(-gamma * (t - tinf))
  
  Ct <- c1 * ((ca * c1a1 * c1a2) + (cb * c1b1 * c1b2) + (cc * c1c1 * c1c2))
  
  c2a1 <- 1 - exp(-alpha * t[t <= tinf])
  c2b1 <- 1 - exp(-beta * t[t <= tinf])
  c2c1 <- 1 - exp(-gamma * t[t <= tinf])
  
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
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.81 p. 48
  Ct <- (dose / dur) * ((A/alpha) * (1 - exp(-alpha * dur)) * exp(-alpha * (t - dur)) +
                          (B/beta) * (1 - exp(-beta * dur)) * exp(-beta * (t - dur)) +
                          (C/gamma) * (1 - exp(-gamma * dur)) * exp(-gamma * (t - dur)))
  Ct[t < dur] <- (dose / dur) * ((A/alpha) * (1 - exp(-alpha * t[t < dur])) +
                                   (B/beta) * (1 - exp(-beta * t[t < dur])) +
                                   (C/gamma) * (1 - exp(-gamma * t[t < dur])))
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
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.84 p. 49
  Ct <- (dose / dur) * ((A/alpha) * (1 - exp(-alpha * dur)) * exp(-alpha * (t - tlag - dur)) +
                          (B/beta) *  (1 - exp(-beta * dur))  * exp(-beta  * (t - tlag - dur)) +
                          (C/gamma) * (1 - exp(-gamma * dur)) * exp(-gamma * (t - tlag - dur)))
  Ct[t <= dur + tlag & t > tlag] <- (dose / dur) * ((A/alpha) * (1 - exp(-alpha * (t[t <= dur + tlag & t > tlag] - tlag))) +
                                                      (B/beta)  * (1 - exp(-beta  * (t[t <= dur + tlag & t > tlag] - tlag))) +
                                                      (C/gamma) * (1 - exp(-gamma * (t[t <= dur + tlag & t > tlag] - tlag))))
  Ct[t <= tlag] <- 0
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
  A <- (1/V1) * (ka / (ka - alpha)) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * (ka / (ka - beta))  * ((k21 - beta)/(beta - alpha))  * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * (ka / (ka - gamma)) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.71 p. 44
  Ct <- dose * ((A * exp(-alpha*t)) + (B * exp(-beta*t)) + (C * exp(-gamma*t)) - ((A + B + C)*exp(-ka*t)))
  Ct
}

#' Calculate C(t) for a 1-compartment linear model
#'
#' @inheritParams calc_sd_1cmt
#' @param ... Passed to `calc_derived_3cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_ss_3cmt <- function(t, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

