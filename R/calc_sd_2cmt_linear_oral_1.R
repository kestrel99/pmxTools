#' Calculate C(t) for a 2-compartment linear model after a single first-order oral dose
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param ka First-order absorption rate constant (/h)
#' @param dose Steady state dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_sd_2cmt_linear_oral_1(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1)
#'
#' @export

calc_sd_2cmt_linear_oral_1 <- function(t, CL, V1, V2, Q, ka, dose) {

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

  ### C(t) after single dose - eq 1.41 p. 25

  Ct <- dose * ((A * exp(-alpha * t)) + (B * exp(-beta * t)) - ((A + B) * exp(-ka * t)))

  Ct
}

