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
  o
}
