#' Calculate derived pharmacokinetic parameters for a 2-compartment linear model
#'
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 Peripheral volume of distribution (L)
#' @param Q Intercompartmental clearance (L/h)
#' @param type Type of half-life to return (\code{"Vss"}, \code{"k10"}, \code{"k12"}, \code{"k21"}, \code{"thalf_alpha"}, \code{"thalf_beta"}, \code{"alpha"}, \code{"beta"}, \code{"trueA"}, \code{"trueB"}, \code{"fracA"}, \code{"fracB"}, \code{"all"}). Default is \code{"all"}).
#' @param sigdig Number of significant digits to be returned. Default is 5.
#' 
#' @return Return a list of derived PK parameters for a 2-compartment linear model given provided clearances and volumes.
#' \itemize{ 
#'   \item \code{Vss}: Volume of distribution at steady state, \eqn{V_{ss}} (L)
#'   \item \code{k10}: First-order elimination rate, \eqn{k_{10}} (/h)
#'   \item \code{k12}: First-order rate of transfer from central to peripheral compartment, \eqn{k_{12}} (/h)
#'   \item \code{k21}: First-order rate of transfer from peripheral to central compartment, \eqn{k_{21}} (/h)   
#'   \item \code{thalf_alpha}: Distributional half-life, \eqn{t_{1/2,\alpha}} (h)
#'   \item \code{thalf_beta}: Terminal half-life, \eqn{t_{1/2,\beta}} (h)
#'   \item \code{alpha}: \eqn{\alpha}
#'   \item \code{beta}: \eqn{\beta}
#'   \item \code{trueA}: true A
#'   \item \code{trueB}: true B
#'   \item \code{fracA}: fractional A
#'   \item \code{fracB}: fractional B
#'  }
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Shafer S. L. \code{CONVERT.XLS}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#'
#' @examples
#' params <- calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
#'
#' @export

calc_derived_2cpt <- function(CL, V1, V2, Q, type="all", sigdig=5) {

  out <- list()
  
  ### variables
  k10   <- CL/V1
  k12   <- Q/V1
  k21   <- Q/V2
  
  out$k10 <- signif(k10, sigdig)
  out$k12 <- signif(k12, sigdig)
  out$k21 <- signif(k21, sigdig)
  
  a0 <- k10 * k21
  a1 <- -(k10 + k12 + k21)
  
  i1 <- (-a1 + sqrt(a1*a1-4*a0))/2
  i2 <- (-a1 - sqrt(a1*a1-4*a0))/2
  
  c1 <- (k21 - i1)/(i2 - i1)/V1
  c2 <- (k21 - i2)/(i1 - i2)/V1

  out$Vss <- signif(V1 + V2, sigdig)
  
  out$thalf_alpha <- signif(log(2)/i1, sigdig)
  out$thalf_beta <- signif(log(2)/i2, sigdig)
  
  out$alpha <- signif(i1, sigdig)
  out$beta  <- signif(i2, sigdig)
  
  out$trueA <- signif(c1, sigdig)
  out$trueB <- signif(c2, sigdig)
  
  out$fracA <- signif(c1*V1, sigdig)
  out$fracB <- signif(c2*V1, sigdig)
  
  if(type=="all") o <- out
  if(type=="k10") o <- out$k10
  if(type=="k12") o <- out$k12
  if(type=="k21") o <- out$k21
  if(type=="trueA") o <- out$trueA
  if(type=="trueB") o <- out$trueB
  if(type=="fracA") o <- out$fracA
  if(type=="fracB") o <- out$fracB
  if(type=="thalf_alpha") o <- out$thalf_alpha
  if(type=="thalf_beta") o <- out$thalf_beta
  if(type=="alpha") o <- out$alpha
  if(type=="beta") o <- out$beta
  if(type=="Vss") o <- out$Vss
  
  o
}

