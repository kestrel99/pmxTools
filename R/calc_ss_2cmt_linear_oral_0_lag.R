#' Calculate C(t) for a 2-compartment linear model at steady-state with zero-order oral dosing and a lag time
#'
#' @param tad Time after last dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param dur Duration of zero-order absorption (h)
#' @param dose Steady state dose
#' @param tau Dosing interval (h)
#' @param tlag Lag time (h)
#'
#' @return Concentration of drug at requested time (\code{tad}) at steady-state, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_ss_2cmt_linear_oral_0_lag(tad = 23, CL = 2.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tau = 24, tlag=2)
#'
#' @export

calc_ss_2cmt_linear_oral_0_lag <- function(tad, CL, V1, V2, Q, dur, dose, tau, tlag) {

  ### microconstants - 1.2 p. 18
  k   <- CL/V1
  k12 <- Q/V1
  k21 <- Q/V2

  ### beta
  beta  <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - (4 * k21 * k)))

  ### alpha
  alpha <- (k21 * k)/beta

  ### macroconstants - 1.2.4 p. 28

  A <- (1/V1) * ((alpha - k21) / (alpha - beta))
  B <- (1/V1) * ((beta - k21) / (beta - alpha))

  ### C(t) at steady state - eq 1.56 p. 33

  c1  <- dose / dur
  ca1 <- A / alpha
  cb1 <- B / beta

  # rest
  
  c1a1 <- 1 - exp(-alpha * dur)
  c1a2 <- exp(-alpha * (tad - tlag - dur))
  c1a3 <- 1 - exp(-alpha * tau)
  
  c1b1 <- 1 - exp(-beta * dur) 
  c1b2 <- exp(-beta * (tad - tlag - dur))
  c1b3 <- 1 - exp(-beta * tau)
  
  Ct <- c1 * (ca1 * ((c1a1 * c1a2) / c1a3) +
              cb1 * ((c1b1 * c1b2) / c1b3))

  # during lag
  
  c2a1 <- 1 - exp(-alpha * dur)
  c2a2 <- exp(-alpha * (tad[tad <= tlag] + tau - tlag - dur))
  c2a3 <- 1 - exp(-alpha * tau)
  
  c2b1 <- 1 - exp(-beta * dur)
  c2b2 <- exp(-beta * (tad[tad <= tlag] + tau - tlag - dur))
  c2b3 <- 1 - exp(-beta * tau)
  
  Ct[tad <= tlag] <- c1 * (ca1 * ((c2a1 * c2a2) / c2a3) +
                           cb1 * ((c2b1 * c2b2) / c2b3))


  # during zero order absorption
  
  c3a1 <- 1 - exp(-alpha * (tad[tad > tlag & tad <= tlag + dur] - tlag))
  c3a2 <- 1 - exp(-alpha * dur)
  c3a3 <- exp(-alpha * (tad[tad > tlag & tad <= tlag + dur] - tlag - dur))
  c3a4 <- 1 - exp(-alpha * tau)
  
  c3b1 <- 1 - exp(-beta * (tad[tad > tlag & tad <= tlag + dur] - tlag))
  c3b2 <- 1 - exp(-beta * dur)
  c3b3 <- exp(-beta * (tad[tad > tlag & tad <= tlag + dur] - tlag - dur))
  c3b4 <- 1 - exp(-beta * tau)
  
  Ct[tad > tlag & tad <= tlag + dur] <-
    c1 * (ca1 * (c3a1 + exp(-alpha * tau) * (c3a2 * c3a3 / c3a4)) +
          cb1 * (c3b1 + exp(-beta * tau) * (c3b2 * c3b3 / c3b4)))

  Ct
}

