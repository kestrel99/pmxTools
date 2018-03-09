#' Calculate C(t) for a 3-compartment linear model at steady-state with first-order oral dosing
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
#'
#' @return Concentration of drug at requested time (\code{t}) at steady state, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_ss_3cmt_linear_oral_1(tad = 11.75, CL = 3.5, V1 = 20,
#'     V2 = 500, V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau = 24)
#'
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

