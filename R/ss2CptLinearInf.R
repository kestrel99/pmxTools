#' Calculate C(t) for a 2-compartment linear model with IV infusion at steady state
#'
#' @param tad Time after last dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param Dose Steady state dose
#' @param Tinf Duration of infusion (h)
#' @param tau Dosing interval (h)
#'
#' @return Concentration of drug at requested time after last dose (\code{tad}), given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://www.lixoft.eu/wp-content/uploads/2015/06/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- ss2CptLinearInf(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, Dose = 10, Tinf = 1, tau = 12)

ss2CptLinearInf <- function(tad, CL, V1, V2, Q, Dose, Tinf, tau) {

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

  ### C(t) at steady state - eq 1.38 p. 23

  Ct <- (Dose/Tinf) * (((A/alpha) * ((((1 - exp(-alpha * Tinf)) * (exp(-alpha * (tad-Tinf)))/(1 - exp(-alpha * tau)))))) +
                         ((B/beta) * ((((1 - exp(-beta * Tinf)) * (exp(-beta * (tad-Tinf)))/(1 - exp(-beta * tau)))))))

  Ct[tad <= Tinf] <- (Dose/Tinf) * (((A/alpha) * (1 - exp(-alpha * tad[tad <= Tinf])) + exp(-alpha * tau) * ((1 - exp(-alpha * Tinf))*(exp(-alpha * (tad[tad <= Tinf] - Tinf)))/(1 - exp(-alpha * tau)))) + ((B/beta) * (1 - exp(-beta * tad[tad <= Tinf])) + exp(-beta * tau) * ((1 - exp(-beta * Tinf))*(exp(-beta * (tad[tad <= Tinf] - Tinf)))/(1 - exp(-beta * tau)))))

  Ct

}

