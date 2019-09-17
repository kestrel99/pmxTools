#' @rdname calc_derived_3cpt
#' @examples
#' params <- calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
#' @export
calc_derived_2cpt <- function(CL, V1, V2, Q2, V=NULL, Q=NULL, type="all", sigdig=5) {
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
  k12   <- Q/V1
  k21   <- Q/V2

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
      fracB=signif(c2*V1, sigdig)
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

