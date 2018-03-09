#' Calculate C(t) for a 2-compartment linear model after a single infusion
#'
#' @param t Time after dose (h)
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param dose Steady state dose
#' @param tinf Duration of infusion (h)
#'
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#'
#' @examples
#' Ctrough <- calc_sd_2cmt_linear_infusion(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10, tinf = 1)
#'
#' @export

calc_sd_2cmt_linear_infusion <- function(t, CL, V1, V2, Q, dose, tinf) {

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

  Ct <- (dose/tinf) * (((A/alpha) * (((1 - exp(-alpha * tinf)) * (exp(-alpha * (t-tinf)))))) +
                       ((B/beta) * (((1 - exp(-beta * tinf)) * (exp(-beta * (t-tinf)))))))

  Ct[t <= tinf] <- (dose/tinf) * (((A/alpha) * (1 - exp(-alpha * t[t <= tinf]))) +
                                    ((B/beta) * (1 - exp(-beta * t[t <= tinf]))))

  Ct
}

