#' Calculate C(t) for a 3-compartment linear model at steady state, with zero-order absorption and lag time
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 First peripheral volume of distribution (L)
#' @param V3 Second peripheral volume of distribution (L)
#' @param Q2 Intercompartmental clearance between V1 and V2 (L/h)
#' @param Q3 Intercompartmental clearance between V2 and V3 (L/h)
#' @param dur Duration of zero-order absorption (h)
#' @param dose Dose
#' @param tau Dosing interval (h)
#' @param tlag Lag time (h)
#'
#' @return Concentration of drug at requested time after dose (\code{tad}) at steady state, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_ss_3cmt_linear_oral_0_lag(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24, tlag = 1.5)
#'
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

