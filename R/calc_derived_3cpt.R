#' Calculate derived pharmacokinetic parameters for a 1-, 2-, or 3-compartment linear model.
#'
#' @param CL Clearance (L/h)
#' @param V1,V Central volume of distribution (L).  Values are synonyms; use only one.
#' @param V2 First peripheral volume of distribution (L)
#' @param V3 Second peripheral volume of distribution (L)
#' @param Q2,Q Intercompartmental clearance from central to first peripheral compartment (L/h).  Values are synonyms; use only one.
#' @param Q3 Intercompartmental clearance from central to second peripheral compartment (L/h)
#' @param type Type of half-life to return. Default is \code{"all"}; see details for other options that are specific to the number of compartments.
#' @param sigdig Number of significant digits to be returned. Default is \code{5}.
#' 
#' @return Return a list of derived PK parameters for a 3-compartment linear model given provided clearances and volumes based on the `type`.
#' \itemize{ 
#'   \item \code{Vss}: Volume of distribution at steady state, \eqn{V_{ss}} (L); 1-, 2-, and 3-compartment
#'   \item \code{k10}: First-order elimination rate, \eqn{k_{10}} (/h); 1-, 2-, and 3-compartment
#'   \item \code{k12}: First-order rate of transfer from central to first peripheral compartment, \eqn{k_{12}} (/h); 2- and 3-compartment
#'   \item \code{k21}: First-order rate of transfer from first peripheral to central compartment, \eqn{k_{21}} (/h); 2- and 3-compartment
#'   \item \code{k13}: First-order rate of transfer from central to second peripheral compartment, \eqn{k_{13}} (/h); 3-compartment
#'   \item \code{k31}: First-order rate of transfer from second peripheral to central compartment,\eqn{k_{31}} (/h); 3-compartment
#'   \item \code{thalf_alpha}: \eqn{t_{1/2,\alpha}} (h); 2- and 3-compartment
#'   \item \code{thalf_beta}: \eqn{t_{1/2,\beta}} (h); 2- and 3-compartment
#'   \item \code{thalf_gamma}: \eqn{t_{1/2,\gamma}} (h); 3-compartment
#'   \item \code{alpha}: \eqn{\alpha}; 1-, 2-, and 3-compartment
#'   \item \code{beta}: \eqn{\beta}; 2- and 3-compartment
#'   \item \code{gamma}: \eqn{\beta}; 3-compartment
#'   \item \code{trueA}: true A; 1-, 2-, and 3-compartment
#'   \item \code{trueB}: true B; 2- and 3-compartment
#'   \item \code{trueC}: true C; 3-compartment
#'   \item \code{fracA}: fractional A; 1-, 2-, and 3-compartment
#'   \item \code{fracB}: fractional B; 2- and 3-compartment
#'   \item \code{fracC}: fractional C; 3-compartment
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

calc_derived_3cpt <- function(CL, V1, V2, V3, Q2, Q3, V=NULL, Q=NULL, type="all", sigdig=5) {
  if (!xor(is.null(V), is.null(V1))) {
    stop("Exactly one of V or V1 may be provided since they are considered synonyms.")
  } else if (!is.null(V)) {
    # Use the preferred synonym
    V1 <- V
  }
  if (!xor(is.null(Q), is.null(Q2))) {
    stop("Exactly one of Q or Q2 may be provided since they are considered synonyms.")
  } else if (!is.null(Q)) {
    # Use the preferred synonym
    Q2 <- Q
  }

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

  out <-
    list(
      k10=signif(k10, sigdig),
      k12=signif(k12, sigdig),
      k21=signif(k21, sigdig),
      k13=signif(k13, sigdig),
      k31=signif(k31, sigdig),
      
      Vss=signif(V1 + V2 + V3, sigdig),
      
      thalf_alpha=signif(log(2)/i1, sigdig),
      thalf_beta =signif(log(2)/i2, sigdig),
      thalf_gamma=signif(log(2)/i3, sigdig),
      
      alpha=signif(i1, sigdig),
      beta =signif(i2, sigdig),
      gamma=signif(i3, sigdig),
      
      trueA=signif(c1, sigdig),
      trueB=signif(c2, sigdig),
      trueC=signif(c3, sigdig),
      
      fracA=signif(c1*V1, sigdig),
      fracB=signif(c2*V1, sigdig),
      fracC=signif(c3*V1, sigdig)
    )
  if(type=="all") {
    o <- out
  } else {
    o <- out[[type]]
    if (is.null(o)) {
      stop("Invalid value for `type`: ", type)
    }
  }
  o
}
