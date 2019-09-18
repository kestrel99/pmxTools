#' Calculate C(t) for a 3-compartment linear model at steady-state
#'
#' @param tad Time after dose (h)
#' @param tau Dosing interval (h)
#' @param dose Dose
#' @param dur Duration of zero-order absorption (h)
#' @param tinf Duration of infusion (h)
#' @param ... Passed to `calc_derived_3cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
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
calc_ss_3cmt_linear_bolus <- function(tad, CL, V1, V2, V3, Q2, Q3, dose, tau) {

  ### microconstants - 1.3 p. 37
  k   <- CL/V1
  k12 <- Q2/V1
  k21 <- Q2/V2
  k13 <- Q3/V1
  k31 <- Q3/V3

  ### a0, a1, a2
  a0 <- k * k21 * k31
  a1 <- (k * k31) + (k21 * k31) + (k21 * k13) + (k * k21) + (k31 * k12)
  a2 <- k + k12 + k13 + k21 + k31

  p  <- a1 - (a2^2)/3
  q  <- ((2 * (a2^3))/27) - ((a1 * a2)/3) + a0
  r1 <- sqrt(-((p^3)/27))
  r2 <- 2 * (r1^(1/3))

  ### phi
  phi   <- acos(-(q/(2 * r1)))/3

  ### alpha
  alpha <- -(cos(phi)*r2 - (a2/3))

  ### beta
  beta  <- -(cos(phi + ((2 * pi)/3)) * r2 - (a2/3))

  ### gamma
  gamma <- -(cos(phi + ((4 * pi)/3)) * r2 - (a2/3))

  ### macroconstants - 1.3.1 p. 39
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.61 p. 40
  Ct <- dose * ((A * (exp(-alpha * tad) / (1 - exp(-alpha * tau)))) +
                (B * (exp(-beta * tad) / (1 - exp(-beta * tau)))) +
                (C * (exp(-gamma * tad) / (1 - exp(-gamma * tau)))))
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady-state with first-order oral dosing with a lag time
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 First peripheral volume of distribution (L)
#' @param V3 Second peripheral volume of distribution (L)
#' @param Q2 Intercompartmental clearance between V1 and V2 (L/h)
#' @param Q3 Intercompartmental clearance between V2 and V3 (L/h)
#' @param ka First-order absorption rate constant (/h)
#' @param dose Dose
#' @param tau Dosing interval (h)
#' @param tlag Lag time (h)
#'
#' @return Concentration of drug at requested time (\code{t}) at steady state after oral dosing, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_ss_3cmt_linear_oral_1_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau=24, tlag = 1.5)
#'
#' @export

calc_ss_3cmt_linear_oral_1_lag <- function(tad, CL, V1, V2, V3, Q2, Q3, ka, dose, tau, tlag) {
  
  ### microconstants - 1.3 p. 37
  k   <- CL/V1
  k12 <- Q2/V1
  k21 <- Q2/V2
  k13 <- Q3/V1
  k31 <- Q3/V3
  
  ### a0, a1, a2
  a0 <- k * k21 * k31
  a1 <- (k * k31) + (k21 * k31) + (k21 * k13) + (k * k21) + (k31 * k12)
  a2 <- k + k12 + k13 + k21 + k31
  
  p  <- a1 - (a2^2)/3
  q  <- ((2 * (a2^3))/27) - ((a1 * a2)/3) + a0
  r1 <- sqrt(-((p^3)/27))
  r2 <- 2 * (r1^(1/3))
  
  ### phi
  phi   <- acos(-(q/(2 * r1)))/3
  
  ### alpha
  alpha <- -(cos(phi)*r2 - (a2/3))
  
  ### beta
  beta  <- -(cos(phi + ((2 * pi)/3)) * r2 - (a2/3))
  
  ### gamma
  gamma <- -(cos(phi + ((4 * pi)/3)) * r2 - (a2/3))
  
  ### macroconstants - 1.3.3 p. 44
  
  A <- (1/V1) * (ka / (ka - alpha)) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * (ka / (ka - beta))  * ((k21 - beta)/(beta - alpha))  * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * (ka / (ka - gamma)) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  
  ### C(t) after single dose with lag time - eq 1.76 p. 45
  
  Ct <- dose * (((A * exp(-alpha*(tad-tlag)))/(1-exp(-alpha*tau))) + ((B * exp(-beta*(tad-tlag)))/(1-exp(-beta*tau))) + ((C * exp(-gamma*(tad-tlag)))/(1-exp(-gamma*tau))) - (((A + B + C)*exp(-ka*(tad-tlag)))/(1-exp(-ka*tau))))
  Ct[tad<tlag] <- dose * (((A * exp(-alpha*(tad[tad<tlag]+tau-tlag)))/(1-exp(-alpha*tau))) + ((B * exp(-beta*(tad[tad<tlag]+tau-tlag)))/(1-exp(-beta*tau))) + ((C * exp(-gamma*(tad[tad<tlag]+tau-tlag)))/(1-exp(-gamma*tau))) - (((A + B + C)*exp(-ka*(tad[tad<tlag]+tau-tlag)))/(1-exp(-ka*tau))))
  
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady state with IV infusions
#' @examples
#' Ct <- calc_ss_3cmt_linear_infusion(tad = 11.75, CL = 2.5, V1 = 20, V2 = 50,
#'     V3 = 100, Q2 = 0.5, Q3 = 0.05, dose = 1000, tinf=1, tau=24)
#' @export
calc_ss_3cmt_linear_infusion <- function(tad, CL, V1, V2, V3, Q2, Q3, dose, tinf, tau) {
  ### microconstants - 1.3 p. 37
  k   <- CL/V1
  k12 <- Q2/V1
  k21 <- Q2/V2
  k13 <- Q3/V1
  k31 <- Q3/V3
  
  ### a0, a1, a2
  a0 <- k * k21 * k31
  a1 <- (k * k31) + (k21 * k31) + (k21 * k13) + (k * k21) + (k31 * k12)
  a2 <- k + k12 + k13 + k21 + k31
  
  p  <- a1 - (a2^2)/3
  q  <- ((2 * (a2^3))/27) - ((a1 * a2)/3) + a0
  r1 <- sqrt(-((p^3)/27))
  r2 <- 2 * (r1^(1/3))
  
  ### phi
  phi   <- acos(-(q/(2 * r1)))/3
  
  ### alpha
  alpha <- -(cos(phi)*r2 - (a2/3))
  
  ### beta
  beta  <- -(cos(phi + ((2 * pi)/3)) * r2 - (a2/3))
  
  ### gamma
  gamma <- -(cos(phi + ((4 * pi)/3)) * r2 - (a2/3))
  
  ### macroconstants - 1.3.2 p. 41
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.68 p. 43
  c1 <- dose / tinf
  ca <- A / alpha
  cb <- B / beta
  cc <- C / gamma
  
  c1a1 <- 1 - exp(-alpha * tinf)
  c1a2 <- exp(-alpha * (tad - tinf))
  c1a3 <- 1 - exp(-alpha * tau)
  
  c1b1 <- 1 - exp(-beta * tinf)
  c1b2 <- exp(-beta * (tad - tinf))
  c1b3 <- 1 - exp(-beta * tau)
  
  c1c1 <- 1 - exp(-gamma * tinf)
  c1c2 <- exp(-gamma * (tad - tinf))
  c1c3 <- 1 - exp(-gamma * tau)
  
  Ct <- c1 * ((ca * ((c1a1 * c1a2)/c1a3)) +
                (cb * ((c1b1 * c1b2)/c1b3)) +
                (cc * ((c1c1 * c1c2)/c1c3)))
  
  c2a1 <- 1 - exp(-alpha * tad[tad <= tinf])
  c2a2 <- exp(-alpha * tau)
  c2a3 <- 1 - exp(-alpha * tinf)
  c2a4 <- exp(-alpha * (tad[tad <= tinf] - tinf))
  c2a5 <- 1 - exp(-alpha * tau)
  
  c2b1 <- 1 - exp(-beta * tad[tad <= tinf])
  c2b2 <- exp(-beta * tau)
  c2b3 <- 1 - exp(-beta * tinf)
  c2b4 <- exp(-beta * (tad[tad <= tinf] - tinf))
  c2b5 <- 1 - exp(-beta * tau)
  
  c2c1 <- 1 - exp(-gamma * tad[tad <= tinf])
  c2c2 <- exp(-gamma * tau)
  c2c3 <- 1 - exp(-gamma * tinf)
  c2c4 <- exp(-gamma * (tad[tad <= tinf] - tinf))
  c2c5 <- 1 - exp(-gamma * tau)
  
  Ct[tad <= tinf] <- c1 * (ca * (c2a1 + (c2a2 * ((c2a3 * c2a4) / c2a5))) +
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
calc_ss_3cmt_linear_oral_0 <- function(tad, CL, V1, V2, V3, Q2, Q3, dur, dose, tau) {
  
  ### microconstants - 1.3 p. 37
  k   <- CL/V1
  k12 <- Q2/V1
  k21 <- Q2/V2
  k13 <- Q3/V1
  k31 <- Q3/V3
  
  ### a0, a1, a2
  a0 <- k * k21 * k31
  a1 <- (k * k31) + (k21 * k31) + (k21 * k13) + (k * k21) + (k31 * k12)
  a2 <- k + k12 + k13 + k21 + k31
  
  p  <- a1 - (a2^2)/3
  q  <- ((2 * (a2^3))/27) - ((a1 * a2)/3) + a0
  r1 <- sqrt(-((p^3)/27))
  r2 <- 2 * (r1^(1/3))
  
  ### phi
  phi   <- acos(-(q/(2 * r1)))/3
  
  ### alpha
  alpha <- -(cos(phi)*r2 - (a2/3))
  
  ### beta
  beta  <- -(cos(phi + ((2 * pi)/3)) * r2 - (a2/3))
  
  ### gamma
  gamma <- -(cos(phi + ((4 * pi)/3)) * r2 - (a2/3))
  
  ### macroconstants - 1.3.4 p. 47
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) after single dose - eq 1.83 p. 49
  Ct <- (dose / dur) * ((A/alpha) * ((1 - exp(-alpha * dur)) * exp(-alpha * (tad - dur)))/(1 - exp(-alpha * tau)) +
                          (B/beta)  * ((1 - exp(-beta * dur)) * exp(-beta * (tad - dur)))/(1 - exp(-beta * tau)) +
                          (C/gamma) * ((1 - exp(-gamma * dur)) * exp(-gamma * (tad - dur)))/(1 - exp(-gamma * tau)))
  
  tm1a <- 1 - exp(-alpha * tad[tad <= dur])
  tm2a <- exp(-alpha * tau)
  tm3a <- 1 - exp(-alpha * dur)
  tm4a <- exp(-alpha * (tad[tad <= dur] - dur))
  tm5a <- 1 - exp(-alpha * tau)
  
  tm1b <- 1 - exp(-beta * tad[tad <= dur])
  tm2b <- exp(-beta * tau)
  tm3b <- 1 - exp(-beta * dur)
  tm4b <- exp(-beta * (tad[tad <= dur] - dur))
  tm5b <- 1 - exp(-beta * tau)
  
  tm1c <- 1 - exp(-gamma * tad[tad <= dur])
  tm2c <- exp(-gamma * tau)
  tm3c <- 1 - exp(-gamma * dur)
  tm4c <- exp(-gamma * (tad[tad <= dur] - dur))
  tm5c <- 1 - exp(-gamma * tau)
  
  Ct[tad <= dur] <- (dose / dur) * ((A/alpha) * (tm1a + (tm2a * ((tm3a * tm4a)/tm5a))) +
                                      (B/beta)  * (tm1b + (tm2b * ((tm3b * tm4b)/tm5b))) +
                                      (C/gamma) * (tm1c + (tm2c * ((tm3c * tm4c)/tm5c))))
  
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady state, with zero-order absorption and lag time
#' @examples
#' Ct <- calc_ss_3cmt_linear_oral_0_lag(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24, tlag = 1.5)
#' @export
calc_ss_3cmt_linear_oral_0_lag <- function(tad, CL, V1, V2, V3, Q2, Q3, dur, dose, tau, tlag) {
  
  ### microconstants - 1.3 p. 37
  k   <- CL/V1
  k12 <- Q2/V1
  k21 <- Q2/V2
  k13 <- Q3/V1
  k31 <- Q3/V3
  
  ### a0, a1, a2
  a0 <- k * k21 * k31
  a1 <- (k * k31) + (k21 * k31) + (k21 * k13) + (k * k21) + (k31 * k12)
  a2 <- k + k12 + k13 + k21 + k31
  
  p  <- a1 - (a2^2)/3
  q  <- ((2 * (a2^3))/27) - ((a1 * a2)/3) + a0
  r1 <- sqrt(-((p^3)/27))
  r2 <- 2 * (r1^(1/3))
  
  ### phi
  phi   <- acos(-(q/(2 * r1)))/3
  
  ### alpha
  alpha <- -(cos(phi)*r2 - (a2/3))
  
  ### beta
  beta  <- -(cos(phi + ((2 * pi)/3)) * r2 - (a2/3))
  
  ### gamma
  gamma <- -(cos(phi + ((4 * pi)/3)) * r2 - (a2/3))
  
  ### macroconstants - 1.3.4 p. 47
  A <- (1/V1) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * ((k21 - beta)/(beta - alpha)) * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) at steady state - eq 1.86 p. 51
  
  tm1a1 <- 1 - exp(-alpha * dur)
  tm2a1 <- exp(-alpha * (tad - tlag - dur))
  tm3a1 <- 1 - exp(-alpha * tau)
  
  tm1b1 <- 1 - exp(-beta * dur)
  tm2b1 <- exp(-beta * (tad - tlag - dur))
  tm3b1 <- 1 - exp(-beta * tau)
  
  tm1c1 <- 1 - exp(-gamma * dur)
  tm2c1 <- exp(-gamma * (tad - tlag - dur))
  tm3c1 <- 1 - exp(-gamma * tau)
  
  Ct <- (dose / dur) * ((A/alpha) * ((tm1a1 * tm2a1)/tm3a1) +
                          (B/beta)  * ((tm1b1 * tm2b1)/tm3b1) +
                          (C/gamma) * ((tm1c1 * tm2c1)/tm3c1))
  
  
  tm1a2 <- 1 - exp(-alpha * (tad[tad > tlag & tad <= tlag + dur] - tlag))
  tm2a2 <- exp(-alpha * tau)
  tm3a2 <- 1 - exp(-alpha * dur)
  tm4a2 <- exp(-alpha * (tad[tad > tlag & tad <= tlag + dur] - tlag - dur))
  tm5a2 <- 1 - exp(-alpha * tau)
  
  tm1b2 <- 1 - exp(-beta * (tad[tad > tlag & tad <= tlag + dur] - tlag))
  tm2b2 <- exp(-beta * tau)
  tm3b2 <- 1 - exp(-beta * dur)
  tm4b2 <- exp(-beta * (tad[tad > tlag & tad <= tlag + dur] - tlag - dur))
  tm5b2 <- 1 - exp(-beta * tau)
  
  tm1c2 <- 1 - exp(-gamma * (tad[tad > tlag & tad <= tlag + dur] - tlag))
  tm2c2 <- exp(-gamma * tau)
  tm3c2 <- 1 - exp(-gamma * dur)
  tm4c2 <- exp(-gamma * (tad[tad > tlag & tad <= tlag + dur] - tlag - dur))
  tm5c2 <- 1 - exp(-gamma * tau)
  
  Ct[tad > tlag & tad <= tlag + dur] <- (dose/dur) * ((A/alpha) * (tm1a2 + (tm2a2 * (tm3a2 * tm4a2)/tm5a2)) +
                                                        (B/beta)  * (tm1b2 + (tm2b2 * (tm3b2 * tm4b2)/tm5b2)) +
                                                        (C/gamma) * (tm1c2 + (tm2c2 * (tm3c2 * tm4c2)/tm5c2)))
  
  tm1a3 <- 1 - exp(-alpha * dur)
  tm2a3 <- exp(-alpha * (tad[tad <= tlag] + tau - tlag - dur))
  tm3a3 <- 1 - exp(-alpha * tau)
  
  tm1b3 <- 1 - exp(-beta * dur)
  tm2b3 <- exp(-beta * (tad[tad <= tlag] + tau - tlag - dur))
  tm3b3 <- 1 - exp(-beta * tau)
  
  tm1c3 <- 1 - exp(-gamma * dur)
  tm2c3 <- exp(-gamma * (tad[tad <= tlag] + tau - tlag - dur))
  tm3c3 <- 1 - exp(-gamma * tau)
  
  Ct[tad <= tlag] <- (dose/dur) * ((A/alpha) * ((tm1a3 * tm2a3) / tm3a3) +
                                     (B/beta)  * ((tm1b3 * tm2b3) / tm3b3) +
                                     (C/gamma) * ((tm1c3 * tm2c3) / tm3c3))
  
  Ct
}

#' @describeIn calc_ss_3cmt
#' Calculate C(t) for a 3-compartment linear model at steady-state with first-order oral dosing
#' @examples
#' Ct <- calc_ss_3cmt_linear_oral_1(tad = 11.75, CL = 3.5, V1 = 20,
#'     V2 = 500, V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau = 24)
#' @export
calc_ss_3cmt_linear_oral_1 <- function(tad, CL, V1, V2, V3, Q2, Q3, ka, dose, tau) {
  
  ### microconstants - 1.3 p. 37
  k   <- CL/V1
  k12 <- Q2/V1
  k21 <- Q2/V2
  k13 <- Q3/V1
  k31 <- Q3/V3
  
  ### a0, a1, a2
  a0 <- k * k21 * k31
  a1 <- (k * k31) + (k21 * k31) + (k21 * k13) + (k * k21) + (k31 * k12)
  a2 <- k + k12 + k13 + k21 + k31
  
  p  <- a1 - (a2^2)/3
  q  <- ((2 * (a2^3))/27) - ((a1 * a2)/3) + a0
  r1 <- sqrt(-((p^3)/27))
  r2 <- 2 * (r1^(1/3))
  
  ### phi
  phi   <- acos(-(q/(2 * r1)))/3
  
  ### alpha
  alpha <- -(cos(phi)*r2 - (a2/3))
  
  ### beta
  beta  <- -(cos(phi + ((2 * pi)/3)) * r2 - (a2/3))
  
  ### gamma
  gamma <- -(cos(phi + ((4 * pi)/3)) * r2 - (a2/3))
  
  ### macroconstants - 1.3.3 p. 44
  A <- (1/V1) * (ka / (ka - alpha)) * ((k21 - alpha)/(alpha - beta)) * ((k31 - alpha)/(alpha - gamma))
  B <- (1/V1) * (ka / (ka - beta))  * ((k21 - beta)/(beta - alpha))  * ((k31 - beta)/(beta - gamma))
  C <- (1/V1) * (ka / (ka - gamma)) * ((k21 - gamma)/(gamma - beta)) * ((k31 - gamma)/(gamma - alpha))
  ### C(t) at steady state - eq 1.73 p. 45
  Ct <- dose * ((A * exp(-alpha*tad))/(1-exp(-alpha*tau)) + (B * exp(-beta*tad))/(1-exp(-beta*tau)) + (C * exp(-gamma*tad))/(1-exp(-gamma*tau)) - ((A + B + C)*exp(-ka*tad))/(1-exp(-ka*tau)))
  Ct
}
