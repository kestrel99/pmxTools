#' Calculate C(t) for a 2-compartment linear model with IV infusion after a single dose
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param Dose Steady state dose
#' @param Tinf Duration of infusion (h)
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- sd2CptLinearInf(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, Dose = 10, Tinf = 1, tau = 12)
#'

sd2CptLinearInf <- function(t, CL, V1, V2, Q, Dose, Tinf, tau) {

  ### microconstants - 1.2 p. 18
  k   <- CL/V1
  k12 <- Q/V1
  k21 <- Q/V2

  ### beta
  beta  <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - (4 * k21 * k)))

  ### alpha
  alpha <- (k21 * k)/beta

  ### macroconstants - 1.2.2 p. 21

  A <- (1/V1) * ((alpha - k21)/(alpha - beta))
  B <- (1/V1) * ((beta - k21)/(beta - alpha))

  ### C(t) after single dose - eq 1.36 p. 22

  Ct <- (Dose/Tinf) * (((A/alpha) * (((1 - exp(-alpha * Tinf)) * (exp(-alpha * (t-Tinf)))))) +
                       ((B/beta) * (((1 - exp(-beta * Tinf)) * (exp(-beta * (t-Tinf)))))))

  Ct[t <= Tinf] <- (Dose/Tinf) * (((A/alpha) * (1 - exp(-alpha * t[t <= Tinf]))) +
                                    ((B/beta) * (1 - exp(-beta * t[t <= Tinf]))))

  Ct
}

