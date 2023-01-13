#' Calculate derived pharmacokinetic parameters for a 1-, 2-, or 3-compartment
#' linear model.
#'
#' @param CL Clearance (volume per time units, e.g. L/h)
#' @param V1,V Central volume of distribution (volume units, e.g. L).  Values
#'   are synonyms; use only one.
#' @param V2 First peripheral volume of distribution (volume units, e.g. L)
#' @param V3 Second peripheral volume of distribution (volume units, e.g. L)
#' @param Q2,Q Intercompartmental clearance from central to first peripheral
#'   compartment (volume per time units, e.g. L/h).  Values are synonyms; use
#'   only one.
#' @param Q3 Intercompartmental clearance from central to second peripheral
#'   compartment (volume per time units, e.g. L/h)
#' @param ka Absorption rate (inverse time units, e.g. 1/h)
#' @param tlag Absorption lag time (time units, e.g. h)
#' @param type Parameters to return. Default is \code{"all"}.  If not
#'   \code{"all"}, this may be a vector from the names of the return value list.
#' @param sigdig Number of significant digits to be returned. Default is
#'   \code{5}.
#' @param verbose For \code{calc_derived()}, provide a message indicating the
#'   type of model detected.
#' @param ... Passed to the other \code{calc_derived_*()} functions.
#' 
#' @return Return a list of derived PK parameters for a 1-, 2-, or 3-compartment
#'   linear model given provided clearances and volumes based on the
#'   \code{type}.
#' \itemize{ 
#'   \item \code{Vss}: Volume of distribution at steady state, \eqn{V_{ss}} (volume units, e.g. L); 1-, 2-, and 3-compartment
#'   \item \code{k10}: First-order elimination rate, \eqn{k_{10}} (inverse time units, e.g. 1/h); 1-, 2-, and 3-compartment
#'   \item \code{k12}: First-order rate of transfer from central to first peripheral compartment, \eqn{k_{12}} (inverse time units, e.g. 1/h); 2- and 3-compartment
#'   \item \code{k21}: First-order rate of transfer from first peripheral to central compartment, \eqn{k_{21}} (inverse time units, e.g. 1/h); 2- and 3-compartment
#'   \item \code{k13}: First-order rate of transfer from central to second peripheral compartment, \eqn{k_{13}} (inverse time units, e.g. 1/h); 3-compartment
#'   \item \code{k31}: First-order rate of transfer from second peripheral to central compartment,\eqn{k_{31}} (inverse time units, e.g. 1/h); 3-compartment
#'   \item \code{thalf_alpha}: \eqn{t_{1/2,\alpha}} (time units, e.g. h); 1-, 2-, and 3-compartment
#'   \item \code{thalf_beta}: \eqn{t_{1/2,\beta}} (time units, e.g. h); 2- and 3-compartment
#'   \item \code{thalf_gamma}: \eqn{t_{1/2,\gamma}} (time units, e.g. h); 3-compartment
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
#' The input parameters with standardized names (\code{V1}, \code{V2},
#' \code{V3}, \code{CL}, \code{Q2}, and \code{Q3}) are also returned in the
#' list, and if provided, additional PK parameters of `ka` and `lag` are also
#' returned in the list.  All inputs may be scalars or vectors.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Shafer S. L. \code{CONVERT.XLS}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and
#'   Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams &
#'   Wilkins, Philadelphia, 2010.
#'
#' @examples
#' params <- calc_derived(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
#' @export
calc_derived <- function(..., verbose=FALSE) {
  arg_names <- names(list(...))
  if (all(arg_names %in% names(formals(calc_derived_1cpt)))) {
    if (verbose) message("Detected 1-compartment model")
    calc_derived_1cpt(...)
  } else if (all(arg_names %in% names(formals(calc_derived_2cpt)))) {
    if (verbose) message("Detected 2-compartment model")
    calc_derived_2cpt(...)
  } else if (all(arg_names %in% names(formals(calc_derived_3cpt)))) {
    if (verbose) message("Detected 3-compartment model")
    calc_derived_3cpt(...)
  } else {
    unknown_params <- setdiff(arg_names, names(formals(calc_derived_3cpt)))
    stop(
      "Could not determine model type based on argument names.  Please check the following argument names: ",
      paste(unknown_params, collapse=", ")
    )
  }
}

#' @rdname calc_derived
#' @examples
#' params <- calc_derived_1cpt(CL=16, V=25)
#' @export
calc_derived_1cpt <- function(CL, V=NULL, V1=NULL, ka=NULL, tlag=NULL, type="all", sigdig=5) {
  if (!xor(is.null(V), is.null(V1))) {
    stop("Exactly one of V or V1 may be provided since they are considered synonyms.")
  } else if (!is.null(V)) {
    # Use the preferred synonym
    V1 <- V
  }
  k10 <- CL/V1
  Vss <- V1
  trueA <- 1/V1
  thalf <- log(2)/k10
  alpha <- k10
  out <-
    list(
      k10=signif(CL/V1, sigdig),
      Vss=signif(Vss, sigdig),
      thalf=signif(thalf, sigdig),
      alpha=signif(alpha, sigdig),
      trueA=signif(1/V1, sigdig),
      fracA=1,
      # Include the macro parameters
      V1=V1,
      CL=CL,
      ka=ka,
      tlag=tlag
    )
  
  if(type=="all") {
    o <- out
  } else {
    o <- out[[type]]
    if (is.null(o)) {
      stop("Invalid value for `type`: ", type)
    }
  }
  o[sapply(o, is.null)] <- NULL   # drops NULL values
  
  if(!is.null(ka)) {
    if(ka < k10) {
      warning("Flip-flop kinetics detected. ka is lower than k10, and half-life may therefore not be similar to the half-life calculated by typical NCA.\n")
    }
  }
  o
}

#' @rdname calc_derived
#' @examples
#' params <- calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
#' @export
calc_derived_2cpt <- function(CL, V1=NULL, V2, Q2=NULL, V=NULL, Q=NULL, ka=NULL, tlag=NULL, type="all", sigdig=5) {
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
  
  a0 <- k10 * k21
  a1 <- -(k10 + k12 + k21)
  
  i1 <- (-a1 + sqrt(a1*a1-4*a0))/2
  i2 <- (-a1 - sqrt(a1*a1-4*a0))/2
  
  c1 <- (k21 - i1)/(i2 - i1)/V1
  c2 <- (k21 - i2)/(i1 - i2)/V1
  
  out <-
    list(
      k10=signif(k10, sigdig),
      k12=signif(k12, sigdig),
      k21=signif(k21, sigdig),
      Vss=signif(V1 + V2, sigdig),
      thalf_alpha=signif(log(2)/i1, sigdig),
      thalf_beta=signif(log(2)/i2, sigdig),
      alpha=signif(i1, sigdig),
      beta=signif(i2, sigdig),
      trueA=signif(c1, sigdig),
      trueB=signif(c2, sigdig),
      fracA=signif(c1*V1, sigdig),
      fracB=signif(c2*V1, sigdig),
      # Include the macro parameters
      V1=V1,
      V2=V2,
      CL=CL,
      Q2=Q2,
      ka=ka,
      tlag=tlag
    )
  
  if(type=="all") {
    o <- out
  } else {
    o <- out[[type]]
    if (is.null(o)) {
      stop("Invalid value for `type`: ", type)
    }
  }
  o[sapply(o, is.null)] <- NULL   # drops NULL values
  
  if(!is.null(ka)) {
    if(ka < k10) {
      warning("Flip-flop kinetics detected. ka is lower than k10, and half-life may therefore not be similar to the half-life calculated by typical NCA.\n")
    }
  }
  
  o
}

#' @rdname calc_derived
#' @examples
#' params <- calc_derived_3cpt(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
#' @export
calc_derived_3cpt <- function(CL, V1=NULL, V2, V3, Q2=NULL, Q3, V=NULL, Q=NULL, ka=NULL, tlag=NULL, type="all", sigdig=5) {
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
  
  i1 <- pmax(root1, root2, root3)
  # There is no pmedian function
  i2 <- mapply(FUN=function(...) median(c(...)), root1, root2, root3)
  i3 <- pmin(root1, root2, root3)
  
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
      fracC=signif(c3*V1, sigdig),
      # Include the macro parameters
      V1=V1,
      V2=V2,
      V3=V3,
      CL=CL,
      Q2=Q2,
      Q3=Q3,
      ka=ka,
      tlag=tlag
    )
  if(type=="all") {
    o <- out
  } else {
    o <- out[[type]]
    if (is.null(o)) {
      stop("Invalid value for `type`: ", type)
    }
  }
  o[sapply(o, is.null)] <- NULL   # drops NULL values
  
  if(!is.null(ka)) {
    if(ka < k10) {
      warning("Flip-flop kinetics detected. ka is lower than k10, and half-life may therefore not be similar to the half-life calculated by typical NCA.\n")
    }
  }
  
  o
}
