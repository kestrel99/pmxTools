#' Calculate C(t) for a 2-compartment linear model at steady-state
#'
#' @param tad Time after dose (h)
#' @param tau Dosing interval (h)
#' @param dose Dose
#' @param dur Duration of zero-order absorption (h)
#' @param tinf Duration of infusion (h)
#' @param ... Passed to `calc_derived_2cpt()`
#' @return Concentration of drug at requested time (\code{t}) after a single dose, given provided set of parameters and variables.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Bertrand J & Mentre F (2008). Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
#' implemented in the Monolix software. \url{http://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#' @export
calc_ss_2cmt <- function(tad, tau, dose, dur=NULL, tinf=NULL, ...) {
  stop("Not yet implemented")
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model with IV bolus dosing at steady-state
#' @examples
#' Ct <- calc_ss_2cmt_linear_bolus(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10, tau=24)
#' @export
calc_ss_2cmt_linear_bolus <- function(tad, CL, V1, V2, Q, dose, tau) {
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

  Ct <- dose * ((A * exp(-alpha * tad) / (1 - exp(-alpha * tau))) + ((B * exp(-beta * tad) / (1 - exp(-beta * tau)))))

  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model with infusion at steady state
#' @examples
#' Ct <- calc_ss_2cmt_linear_infusion(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 10, tinf = 1, tau = 12)
#' @export
calc_ss_2cmt_linear_infusion <- function(tad, CL, V1, V2, Q, dose, tinf, tau) {
  
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
  InfA2 <- exp(-1*alpha*tau)*(((1-exp(-1*alpha*tinf))* exp(-1*alpha*(Time-tinf)))/(1-exp(-1*alpha*tau)))
  InfB1 <- 1-exp(-1*beta*Time)
  InfB2 <- exp(-1*beta*tau)*(((1-exp(-1*beta*tinf))* exp(-1*beta*(Time-tinf)))/(1-exp(-1*beta*tau)))
  
  PInfA2 <- (((1-exp(-1*alpha*tinf))* exp(-1*alpha*(Time-tinf)))/(1-exp(-1*alpha*tau)))
  PInfB2 <- (((1-exp(-1*beta *tinf))* exp(-1*beta *(Time-tinf)))/(1-exp(-1*beta *tau)))
  
  InfProf  <- A/alpha*(InfA1+InfA2)+B/beta*(InfB1+InfB2)
  PInfProf <- A/alpha*(PInfA2)+B/beta*(PInfB2)
  
  Ct <- dose/tinf*(InfProf*(Time<=tinf)+PInfProf*(Time>tinf))
  
  Ct
  
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with zero-order oral dosing
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_0(tad = 23, CL = 2.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tau = 24)
#' @export
calc_ss_2cmt_linear_oral_0 <- function(tad, CL, V1, V2, Q, dur, dose, tau) {
  
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
  
  ### C(t) at steady state - eq 1.53 p. 30
  
  c1  <- dose / dur
  ca1 <- A / alpha
  cb1 <- B / beta
  
  c1a1 <- 1 - exp(-alpha * dur)
  c1a2 <- exp(-alpha * (tad - dur))
  c1a3 <- 1 - exp(-alpha * tau)
  
  c1b1 <- 1 - exp(-beta * dur)
  c1b2 <- exp(-beta * (tad - dur))
  c1b3 <- 1 - exp(-beta * tau)
  
  Ct <- c1 * ((ca1 * ((c1a1 * c1a2) / c1a3)) +
                (cb1 * ((c1b1 * c1b2) / c1b3)))
  
  
  c2a1 <- 1 - exp(-alpha * tad[tad < dur])
  c2a2 <- exp(-alpha * tau)
  c2a3 <- 1 - exp(-alpha * dur)
  c2a4 <- exp(-alpha * (tad[tad < dur] - dur))
  c2a5 <- 1 - exp(-alpha * tau)
  
  c2b1 <- 1 - exp(-beta * tad[tad < dur])
  c2b2 <- exp(-beta * tau)
  c2b3 <- 1 - exp(-beta * dur)
  c2b4 <- exp(-beta * (tad[tad < dur] - dur))
  c2b5 <- 1 - exp(-beta * tau)
  
  Ct[tad < dur] <- c1 * ((ca1 * (c2a1 + c2a2*((c2a3 * c2a4)/c2a5))) +
                           (cb1 * (c2b1 + c2b2*((c2b3 * c2b4)/c2b5))))
  
  Ct
}

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with first-order oral dosing
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_1_lag(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tau=24, tlag=2)
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

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with zero-order oral dosing and a lag time
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_0_lag(tad = 23, CL = 2.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, dur = 1, tau = 24, tlag=2)
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

#' @describeIn calc_ss_2cmt
#' Calculate C(t) for a 2-compartment linear model at steady-state with first-order oral dosing
#' @examples
#' Ct <- calc_ss_2cmt_linear_oral_1(tad = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
#'     dose = 1000, ka = 1, tau=24)
#' @export
calc_ss_2cmt_linear_oral_1 <- function(tad, CL, V1, V2, Q, ka, dose, tau) {
  
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
  ### C(t) after single dose - eq 1.43 p. 25
  Ct <- dose * (((A * exp(-alpha * tad)) / (1 - exp(-alpha * tau)))  + ((B * exp(-beta * tad)) / (1 - exp(-beta * tau))) - (((A + B) * exp(-ka * tad)) / (1 - exp(-ka * tau))))
  Ct
}
