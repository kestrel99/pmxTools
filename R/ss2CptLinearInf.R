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

  Time <- tad
  
  ### microconstants - 1.2 p. 18
  k12 <- Q/V1
  k21 <- Q/V2
  k   <- CL/V1

  ### additional derivations
  agg     <- k12 + k21 + k
  beta    <- 0.5*(agg - sqrt(agg*agg - 4*k21*k))
  alpha   <- k21*k/beta
  t_alpha <- log(2)/alpha
  t_beta  <- log(2)/beta

  ### macroconstants - 1.2.2 p. 21
  A <- 1 / V1 * ((alpha-k21)/(alpha-beta))
  B <- 1 / V1 * ((beta-k21)/(beta-alpha))

  ### C(t) at steady state - eq 1.38 p. 23

  InfA1 <- 1-exp(-1*alpha*Time)
  InfA2 <- exp(-1*alpha*tau)*(((1-exp(-1*alpha*Tinf))* exp(-1*alpha*(Time-Tinf)))/(1-exp(-1*alpha*tau)))
  InfB1 <- 1-exp(-1*beta*Time)
  InfB2 <- exp(-1*beta*tau)*(((1-exp(-1*beta*Tinf))* exp(-1*beta*(Time-Tinf)))/(1-exp(-1*beta*tau)))
  
  PInfA2 <- (((1-exp(-1*alpha*Tinf))* exp(-1*alpha*(Time-Tinf)))/(1-exp(-1*alpha*tau)))
  PInfB2 <- (((1-exp(-1*beta *Tinf))* exp(-1*beta *(Time-Tinf)))/(1-exp(-1*beta *tau)))
  
  InfProf  <- A/alpha*(InfA1+InfA2)+B/beta*(InfB1+InfB2)
  PInfProf <- A/alpha*(PInfA2)+B/beta*(PInfB2)
  
  Ct <- Dose/Tinf*(InfProf*(Time<=Tinf)+PInfProf*(Time>Tinf))
    
  Ct

}

