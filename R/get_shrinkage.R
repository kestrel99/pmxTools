#' Extract shrinkage estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param output A flag specifying the shrinkage estimates to be output. Valid flag
#' values are \code{eta} (the default), \code{epsilon}, or \code{all}.
#' @param type Specifies the type of shrinkage to report. Valid values are \code{sd} 
#' (standard deviation, the default) or \code{vr} (variance, if present in the XML output). 
#' @param sigdig Specifies the number of significant digits to be provided (default=3).
#' @param est.step Specifies which estimation step to return parameters from (default is the last).
#'
#' @return A named vector of NONMEM shrinkage estimates, or in the case of \code{all},
#' a list of named vectors.
#'
#' \code{eta} returns a vector of ETA shrinkages, as reported by NONMEM.
#' \code{epsilon} returns EPSILON shrinkage, as reported by NONMEM.
#' \code{all} returns both ETA and EPSILON shrinkage estimates as a list of vectors.
#' 
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  nmOutput <- read_nm("run315.xml")
#'  shr <- get_shrinkage(nmOutput, output="all")
#' }
#'
#' @export

get_shrinkage <- function (x, output = "eta", type="sd", sigdig = 3, est.step=NULL) {
  if (!(output %in% c("eta", "epsilon", "all"))) {
    stop("Please select a valid output option (eta, epsilon, all).")
  }
  
  if(is.null(est.step)) {
    no_steps <- sum(stringr::str_count(names(x$nonmem$problem), "estimation"))
  } else {
    no_steps <- est.step
  }
  
  ind_est  <- match("estimation", names(x$nonmem$problem))-1+no_steps
  
  if(!is.null(x$nonmem$problem[[ind_est]]$etashrink)) {
    eta <- signif(as.numeric(x$nonmem$problem[[ind_est]]$etashrink[seq(1, 
                                                                       length(x$nonmem$problem[[ind_est]]$etashrink) - 1, by = 2)]), 
                  sigdig)
    names(eta) <- as.character(x$nonmem$problem[[ind_est]]$etashrink[seq(2, 
                                                                         length(x$nonmem$problem[[ind_est]]$etashrink) - 1, by = 2)])
    eps <- signif(as.numeric(x$nonmem$problem[[ind_est]]$epsshrink[seq(1, 
                                                                       length(x$nonmem$problem[[ind_est]]$epsshrink) - 1, by = 2)]), 
                  sigdig)
    names(eps) <- as.character(x$nonmem$problem[[ind_est]]$epsshrink[seq(2, 
                                                                         length(x$nonmem$problem[[ind_est]]$epsshrink) - 1, by = 2)])
  } else {
    eta <- signif(as.numeric(x$nonmem$problem[[ind_est]]$etashrinksd[seq(1, 
                                                                         length(x$nonmem$problem[[ind_est]]$etashrinksd) - 1, by = 2)]), 
                  sigdig)
    names(eta) <- as.character(x$nonmem$problem[[ind_est]]$etashrinksd[seq(2, 
                                                                           length(x$nonmem$problem[[ind_est]]$etashrinksd) - 1, by = 2)])
    eps <- signif(as.numeric(x$nonmem$problem[[ind_est]]$epsshrinksd[seq(1, 
                                                                         length(x$nonmem$problem[[ind_est]]$epsshrinksd) - 1, by = 2)]), 
                  sigdig)
    names(eps) <- as.character(x$nonmem$problem[[ind_est]]$epsshrinksd[seq(2, 
                                                                           length(x$nonmem$problem[[ind_est]]$epsshrinksd) - 1, by = 2)]) 
    etavr <- signif(as.numeric(x$nonmem$problem[[ind_est]]$etashrinkvr[seq(1, 
                                                                           length(x$nonmem$problem[[ind_est]]$etashrinkvr) - 1, by = 2)]), 
                    sigdig)
    names(etavr) <- as.character(x$nonmem$problem[[ind_est]]$etashrinkvr[seq(2, 
                                                                             length(x$nonmem$problem[[ind_est]]$etashrinkvr) - 1, by = 2)])
    epsvr <- signif(as.numeric(x$nonmem$problem[[ind_est]]$epsshrinkvr[seq(1, 
                                                                           length(x$nonmem$problem[[ind_est]]$epsshrinkvr) - 1, by = 2)]), 
                    sigdig)
    names(epsvr) <- as.character(x$nonmem$problem[[ind_est]]$epsshrinkvr[seq(2, 
                                                                             length(x$nonmem$problem[[ind_est]]$epsshrinkvr) - 1, by = 2)]) 
  }
  if (output == "eta" & type=="sd") {
    out <- eta
  }
  if (output == "eta" & type=="vr") {
    out <- etavr
  }
  
  if (output == "epsilon" & type=="sd") {
    out <- eps
  }
  if (output == "epsilon" & type=="vr") {
    out <- epsvr
  }  
  
  if (output == "all") {
    if(!is.null(x$nonmem$problem[[ind_est]]$etashrink)) {
      out <- list()
      out$etasd <- eta
      out$etavr <- etavr
      out$epsilonsd <- eps
      out$epsilonvr <- epsvr
    } else {
      out <- list()
      out$etasd <- eta
      out$epsilonsd <- eps
    }
  } 
  out
}
