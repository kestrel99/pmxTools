#' Calculate C(t) for a 3-compartment linear model at steady state with IV infusions
#'
#' @param tad Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 First peripheral volume of distribution (L)
#' @param V3 Second peripheral volume of distribution (L)
#' @param Q2 Intercompartmental clearance between V1 and V2 (L/h)
#' @param Q3 Intercompartmental clearance between V2 and V3 (L/h)
#' @param dose Dose
#' @param tinf Duration of infusion (h)
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_ss_3cmt_linear_infusion(tad = 11.75, CL = 2.5, V1 = 20, V2 = 50, V3 = 100, Q2 = 0.5, Q3 = 0.05,
#'     dose = 1000, tinf=1, tau=24)
#'
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

