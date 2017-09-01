#' Calculate C(t) for a 3-compartment linear model at steady state, with zero-order absorption
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
#'
#' @return Concentration of drug at requested time after dose (\code{tad}) at steady state, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_ss_3cmt_linear_oral_0(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
#'     V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24)
#'
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

