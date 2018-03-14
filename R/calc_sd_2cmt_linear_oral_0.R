#' Calculate C(t) for a 2-compartment linear model after a single zero-order oral dose
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param dur Duration of zero-order absorption (h)
#' @param dose Steady state dose
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.
#' 
#' @examples
#' Ctrough <- calc_sd_2cmt_linear_oral_0(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1)
#'
#' @export

calc_sd_2cmt_linear_oral_0 <- function(t, CL, V1, V2, Q, dur, dose) {

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

  ### C(t) after single dose - eq 1.51 p. 29

  Ct <- (dose / dur) * (((A / alpha) * (1 - exp(-alpha * dur)) * exp(-alpha * (t - dur))) +
                      ((B / beta) * (1 - exp(-beta * dur)) * exp(-beta * (t - dur))))

  Ct[t < dur] <- (dose / dur) * ((A / alpha) * (1 - exp(-alpha * t[t < dur])) +
                                 (B / beta) * (1 - exp(-beta * t[t < dur])))

  Ct
}

