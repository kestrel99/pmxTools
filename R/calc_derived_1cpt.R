#' Calculate derived pharmacokinetic parameters for a 1-compartment linear model
#'
#' @param CL Clearance (L/h)
#' @param V  Central volume of distribution (L)
#' @param type Type of half-life to return (\code{"Vss"}, \code{"k10"}, \code{"thalf"}, \code{"alpha"}, \code{"trueA"}, \code{"fracA"}, \code{"all"}). Default is \code{"all"}).
#' @param sigdig Number of significant digits to be returned. Default is 5.
#' 
#' @return Return a list of derived PK parameters for a 2-compartment linear model given provided clearances and volumes.
#' \itemize{ 
#'   \item \code{Vss}: \eqn{V_{ss}} (L)
#'   \item \code{k10}: \eqn{k_{10}} (/h)
#'   \item \code{thalf}: \eqn{t_{1/2}} (h)
#'   \item \code{alpha}: \eqn{\alpha}
#'   \item \code{trueA}: true A
#'   \item \code{fracA}: fractional A
#'  }
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Shafer S. L. \code{CONVERT.XLS}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#'
#' @examples
#' params <- calc_derived_1cpt(CL=16, V=25)
#'
#' @export

calc_derived_1cpt <- function(CL, V, type="all", sigdig=5) {

  out <- list()
  
  ### variables
  out$k10   <- signif(CL/V, sigdig)

  out$Vss <- V
  
  out$thalf <- signif(log(2)/out$k10, sigdig)
  
  out$alpha <- out$k10
  
  out$trueA <- signif(1/V, sigdig)
  out$fracA <- 1
  
  
  if(type=="all") o <- out
  if(type=="k10") o <- out$k10
  if(type=="trueA") o <- out$trueA
  if(type=="fracA") o <- out$fracA
  if(type=="thalf") o <- out$thalf
  if(type=="alpha") o <- out$alpha
  if(type=="Vss") o <- out$Vss
  
  o

  
}

