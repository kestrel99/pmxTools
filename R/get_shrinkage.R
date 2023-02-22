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
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
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
  
  eta <- NA
  eps <- NA
  etavr <- NA
  epsvr <- NA
  
  if(is.null(est.step)) {
    no_steps <- sum(stringr::str_count(names(x$nonmem$problem), "estimation"))
  } else {
    no_steps <- est.step
  }
  
  ind_est  <- match("estimation", names(x$nonmem$problem))-1+no_steps
  
  if(!is.null(x$nonmem$problem[[ind_est]]$etashrink)) {
    
    eta <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$etashrink)), sigdig)
    names(eta) <- paste("ETA",1:length(eta),sep = "")
    
    eps <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$epsshrink)), sigdig)
    names(eps) <- paste("EPS",1:length(eps),sep = "")
    
  } else {
    
    eta <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$etashrinksd)), sigdig)
    names(eta) <- paste("ETA",1:length(eta),sep = "")
    
    eps <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$epsshrinksd)), sigdig)
    names(eps) <- paste("EPS",1:length(eps),sep = "")
    
    etavr <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$etashrinkvr)), sigdig)
    names(etavr) <- paste("ETA",1:length(eta),sep = "")
  
    epsvr <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$epsshrinkvr)), sigdig)
    names(epsvr) <- paste("EPS",1:length(eps),sep = "")        
    
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
      out <- list()
      out$etasd <- eta
      out$etavr <- etavr
      out$epsilonsd <- eps
      out$epsilonvr <- epsvr
  } 
  out
}
