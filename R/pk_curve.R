#' Provide concentration-time curves.
#'
#' @param t Observation time in h, specified as a vector.
#' @param model The model to use. Must be one of "1cmt_bolus", "1cmt_infusion", "1cmt_oral", "2cmt_bolus", "2cmt_infusion",
#' "2cmt_oral", "3cmt_bolus", "3cmt_infusion", "3cmt_oral". The default is "1cmt_oral".
#' @param params A named list containing parameter values for the selected model type.
#' @param dose Dose amount.
#' @param ii Interdose interval (or tau), in hours (default 24).
#' @param addl Number of additional doses (default 0).
#' @param ss Assume steady state concentration (default \code{FALSE}).
#'
#' @return A data frame containing times (\code{t}) and concentrations (\code{cp}).
#' @export
#'
#' @examples
#' plot(pk_curve(t=seq(0,72,by=0.1), model="3cmt_oral", ii=12, addl=5,
#'   params=list(CL=2.5, V1=25, V2=2, V3=5, Q2=0.5, Q3=0.25, ka=1)), type="l")

pk_curve <- function(t, model = "1cmt_oral", params=list(ka=2.77, CL=2.5, V=25), dose=600, ii=24, addl=0, ss=F) {

  if(model == "1cmt_oral") {
    if (!all(sort(names(params)) == sort(c("ka","CL","V")))) {
      stop("The one-compartment oral model requires parameter values for CL, V and ka.")
    }
    baseCurve <- calc_sd_1cmt_linear_oral_1(t=t, CL=params$CL, V=params$V, ka=params$ka, dose=dose)
  }

  if(model == "1cmt_bolus") {
    if (!all(sort(names(params)) == sort(c("CL","V")))) {
      stop("The one-compartment IV bolus model requires parameter values for CL and V.")
    }
    baseCurve <- calc_sd_1cmt_linear_bolus(t=t, CL=params$CL, V=params$V, dose=dose)
  }

  if(model == "1cmt_infusion") {
    if (!all(sort(names(params)) == sort(c("CL","V","tinf")))) {
      stop("The one-compartment IV infusion model requires parameter values for CL, V and tinf.")
    }
    baseCurve <- calc_sd_1cmt_linear_infusion(t=t, CL=params$CL, V=params$V, tinf=params$tinf, dose=dose)
  }

  if(model == "2cmt_oral") {
    if (!all(sort(names(params)) == sort(c("ka","CL","V1","V2","Q")))) {
      stop("The two-compartment oral model requires parameter values for CL, V1, V2, Q and ka.")
    }
    baseCurve <- calc_sd_2cmt_linear_oral_1(t=t, CL=params$CL, V1=params$V1, V2=params$V2, Q=params$Q, ka=params$ka, dose=dose)
  }

  if(model == "2cmt_bolus") {
    if (!all(sort(names(params)) == sort(c("CL","V1","V2","Q")))) {
      stop("The one-compartment IV bolus model requires parameter values for CL, V1, V2 and Q.")
    }
    baseCurve <- calc_sd_2cmt_linear_bolus(t=t, CL=params$CL, V1=params$V1, V2=params$V2, Q=params$Q, dose=dose)
  }

  if(model == "2cmt_infusion") {
    if (!all(sort(names(params)) == sort(c("CL","V1","V2","Q","tinf")))) {
      stop("The two-compartment IV infusion model requires parameter values for CL, V1, V2, Q and tinf.")
    }
    baseCurve <- calc_sd_2cmt_linear_infusion(t=t, CL=params$CL, V1=params$V1, V2=params$V2,
                                              Q=params$Q, tinf=params$tinf, dose=dose)
  }
  if(model == "3cmt_oral") {
    if (!all(sort(names(params)) == sort(c("ka","CL","V1","V2","V3","Q2","Q3")))) {
      stop("The three-compartment oral model requires parameter values for CL, V1, V2, V3, Q2, Q3 and ka.")
    }
    baseCurve <- calc_sd_3cmt_linear_oral_1(t=t, CL=params$CL, V1=params$V1, V2=params$V2, V3=params$V3,
                                            Q2=params$Q2, Q3=params$Q3, ka=params$ka, dose=dose)
  }

  if(model == "3cmt_bolus") {
    if (!all(sort(names(params)) == sort(c("CL","V1","V2","V3","Q2","Q3")))) {
      stop("The three-compartment IV bolus model requires parameter values for CL, V1, V2, V3, Q2 and Q3.")
    }
    baseCurve <- calc_sd_3cmt_linear_bolus(t=t, CL=params$CL, V1=params$V1, V2=params$V2, V3=params$V3,
                                           Q2=params$Q2, Q3=params$Q3, dose=dose)
  }

  if(model == "3cmt_infusion") {
    if (!all(sort(names(params)) == sort(c("CL","V1","V2","V3","Q2","Q3","tinf")))) {
      stop("The three-compartment IV infusion model requires parameter values for CL, V1, V2, V3, Q2, Q3 and tinf.")
    }
    baseCurve <- calc_sd_3cmt_linear_infusion(t=t, CL=params$CL, V1=params$V1, V2=params$V2, V3=params$V3,
                                              Q2=params$Q2, Q3=params$Q3, tinf=params$tinf, dose=dose)
  }

  out <- data.frame(t  = t,
                    cp = baseCurve)

  if(addl>0) {
    for (n in 1:addl) {
      zeroes <- rep(0, times=length(t[t<(ii*n)]))
      addCurve <- c(zeroes, baseCurve)
      addCurve <- head(addCurve, length(t))
      out$cp <- out$cp + addCurve
    }
  }

  out

}
