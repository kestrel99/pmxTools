#' Calculate derived pharmacokinetic parameters for a 3-compartment linear model
#'
#' @param CL Clearance (L/h)
#' @param V1 Central volume of distribution (L)
#' @param V2 First peripheral volume of distribution (L)
#' @param V3 Second peripheral volume of distribution (L)
#' @param Q2 Intercompartmental clearance from central to first peripheral compartment (L/h)
#' @param Q3 Intercompartmental clearance from central to second peripheral compartment (L/h)
#' @param type Type of half-life to return (\code{"Vss"}, \code{"k10"}, \code{"k12"}, \code{"k13"}, \code{"k21"}, \code{"k31"}, \code{"thalf_alpha"}, \code{"thalf_beta"}, \code{"thalf_gamma"}, \code{"alpha"}, \code{"beta"}, \code{"gamma"}, \code{"trueA"}, \code{"trueB"}, \code{"trueC"}, \code{"fracA"}, \code{"fracB"}, \code{"fracC"}, \code{"all"}). Default is \code{"all"}).
#' @param sigdig Number of significant digits to be returned. Default is 5.
#' 
#' @return Return a list of derived PK parameters for a 3-compartment linear model given provided clearances and volumes.
#' \itemize{ 
#'   \item \code{Vss}: Volume of distribution at steady state, \eqn{V_{ss}} (L)
#'   \item \code{k10}: First-order elimination rate, \eqn{k_{10}} (/h)
#'   \item \code{k12}: First-order rate of transfer from central to first peripheral compartment, \eqn{k_{12}} (/h)
#'   \item \code{k21}: First-order rate of transfer from first peripheral to central compartment, \eqn{k_{21}} (/h)   
#'   \item \code{k13}: First-order rate of transfer from central to second peripheral compartment, \eqn{k_{13}} (/h)
#'   \item \code{k31}: First-order rate of transfer from second peripheral to central compartment,\eqn{k_{31}} (/h)  
#'   \item \code{thalf_alpha}: \eqn{t_{1/2,\alpha}} (h)
#'   \item \code{thalf_beta}: \eqn{t_{1/2,\beta}} (h)
#'   \item \code{thalf_gamma}: \eqn{t_{1/2,\gamma}} (h)
#'   \item \code{alpha}: \eqn{\alpha}
#'   \item \code{beta}: \eqn{\beta}
#'   \item \code{gamma}: \eqn{\beta}
#'   \item \code{trueA}: true A
#'   \item \code{trueB}: true B
#'   \item \code{trueC}: true C
#'   \item \code{fracA}: fractional A
#'   \item \code{fracB}: fractional B
#'   \item \code{fracC}: fractional C
#'  }
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Shafer S. L. \code{CONVERT.XLS}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#'
#' @examples
#' params <- calc_derived_3cpt(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
#'
#' @export

calc_derived_3cpt <- function(CL, V1, V2, V3, Q2, Q3, type="all", sigdig=5) {

  out <- list()
  
  ### variables
  k10   <- CL/V1
  k12   <- Q2/V1
  k21   <- Q2/V2
  k13   <- Q3/V1
  k31   <- Q3/V3
  
  a0 <- k10 * k21 * k31
  a1 <- (k10 * k31) + (k21 * k31) + (k21 * k13) + (k10 * k21) + (k31 * k12)
  a2 <- k10 + k12 + k13 + k21 + k31
    
  p   <- a1 - (a2 * a2 / 3)
  q   <- (2 * a2 * a2 * a2 / 27) - (a1 * a2 / 3) + a0
  r1  <- sqrt(-(p * p * p)/27)
  phi <- acos((-q/2)/r1)/3
  r2  <- 2 * exp(log(r1)/3)
  
  root1 <- -(cos(phi) * r2 - a2/3)
  root2 <- -(cos(phi + 2 * pi/3) * r2 - a2/3)
  root3 <- -(cos(phi + 4 * pi/3) * r2 - a2/3)
    
  i1 <- max(c(root1, root2, root3))
  i2 <- median(c(root1, root2, root3))
  i3 <- min(c(root1, root2, root3))
  
  c1 <- (k21 - i1) * (k31 - i1) / (i1 - i2) / (i1 - i3) / V1
  c2 <- (k21 - i2) * (k31 - i2) / (i2 - i1) / (i2 - i3) / V1
  c3 <- (k21 - i3) * (k31 - i3) / (i3 - i2) / (i3 - i1) / V1

  out$k10 <- signif(k10, sigdig)
  out$k12 <- signif(k12, sigdig)
  out$k21 <- signif(k21, sigdig)
  out$k13 <- signif(k13, sigdig)
  out$k31 <- signif(k31, sigdig)
  
  out$Vss <- signif(V1 + V2 + V3, sigdig)
  
  out$thalf_alpha <- signif(log(2)/i1, sigdig)
  out$thalf_beta  <- signif(log(2)/i2, sigdig)
  out$thalf_gamma <- signif(log(2)/i3, sigdig)
  
  out$alpha <- signif(i1, sigdig)
  out$beta  <- signif(i2, sigdig)
  out$gamma <- signif(i3, sigdig)
  
  out$trueA <- signif(c1, sigdig)
  out$trueB <- signif(c2, sigdig)
  out$trueC <- signif(c3, sigdig)
  
  out$fracA <- signif(c1*V1, sigdig)
  out$fracB <- signif(c2*V1, sigdig)
  out$fracC <- signif(c3*V1, sigdig)
  
  if(type=="all") o <- out
  if(type=="k10") o <- out$k10
  if(type=="k12") o <- out$k12
  if(type=="k21") o <- out$k21
  if(type=="k13") o <- out$k13
  if(type=="k31") o <- out$k31
  if(type=="trueA") o <- out$trueA
  if(type=="trueB") o <- out$trueB
  if(type=="trueC") o <- out$trueC
  if(type=="fracA") o <- out$fracA
  if(type=="fracB") o <- out$fracB
  if(type=="fracC") o <- out$fracC
  if(type=="thalf_alpha") o <- out$thalf_alpha
  if(type=="thalf_beta") o <- out$thalf_beta
  if(type=="thalf_gamma") o <- out$thalf_gamma
  if(type=="alpha") o <- out$alpha
  if(type=="beta") o <- out$beta
  if(type=="gamma") o <- out$gamma
  if(type=="Vss") o <- out$Vss
  
  o
  
}

