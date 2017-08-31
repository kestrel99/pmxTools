#' Calculate C(t) for a 2-compartment linear model at steady-state with first-order oral dosing
#'
#' @param tad Time after last dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param ka First-order absorption rate constant (/h)
#' @param dose Steady state dose
#' @param tau Dosing interval (h)
#' @param tlag Lag time (h)
#'
#' @return Concentration of drug at requested time after last dose (\code{tad}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_ss_2cmt_linear_oral_1_lag(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tau=24, tlag=2)
#'
#' @export

calc_ss_2cmt_linear_oral_1_lag <- function(tad, CL, V1, V2, Q, ka, dose, tau, tlag) {

  ### microconstants - 1.2 p. 18
  k   <- CL/V1
  k12 <- Q/V1
  k21 <- Q/V2

  ### beta
  beta  <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k)^2 - (4 * k21 * k)))

  ### alpha
  alpha <- (k21 * k)/beta

  ### macroconstants - 1.2.3 p. 24

  A <- (ka/V1) * ((k21 - alpha) / ((ka - alpha) * (beta - alpha)))
  B <- (ka/V1) * ((k21 - beta) / ((ka - beta) * (alpha - beta)))

  ### C(t) after single dose - eq 1.46 p. 26

  Ct <- dose * (((A * exp(-alpha * (tad - tlag))) / (1 - exp(-alpha * tau)))  +
                  ((B * exp(-beta * (tad - tlag))) / (1 - exp(-beta * tau))) -
                  (((A + B) * exp(-ka * (tad - tlag))) / (1 - exp(-ka * tau))))


  Ct[tad < tlag] <- dose * (((A * exp(-alpha * (tad[tad < tlag] + tau - tlag))) / (1 - exp(-alpha * tau)))  +
                              ((B * exp(-beta * (tad[tad < tlag] + tau - tlag))) / (1 - exp(-beta * tau))) -
                              (((A + B) * exp(-ka * (tad[tad < tlag] + tau - tlag))) / (1 - exp(-ka * tau))))

  Ct
}

