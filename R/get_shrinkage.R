#' Extract shrinkage estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param output A flag specifying the shrinkage estimates to be output. Valid flag
#' values are \code{eta} (the default), \code{epsilon}, or \code{all}.
#' @param type Specifies the type of shrinkage to report. Valid values are \code{sd} 
#' (standard deviation, the default) or \code{vr} (variance, if present in the XML output). 
#' @param sigdig Specifies the number of significant digits to be provided (default=3).
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

get_shrinkage <- function (x, output = "eta", type="sd", sigdig = 3) {
  if (!(output %in% c("eta", "epsilon", "all"))) {
    stop("Please select a valid output option (eta, epsilon, all).")
  }
  if (!(type %in% c("sd", "vr"))) {
    stop("Please select a valid shrinkage type option (sd, vr).")
  }
  
  if(!is.null(x$nonmem$problem$estimation$etashrink)) {
    eta <- signif(as.numeric(x$nonmem$problem$estimation$etashrink[seq(1, 
                                                                       length(x$nonmem$problem$estimation$etashrink) - 1, by = 2)]), 
                  sigdig)
    names(eta) <- as.character(x$nonmem$problem$estimation$etashrink[seq(2, 
                                                                         length(x$nonmem$problem$estimation$etashrink) - 1, by = 2)])
    eps <- signif(as.numeric(x$nonmem$problem$estimation$epsshrink[seq(1, 
                                                                       length(x$nonmem$problem$estimation$epsshrink) - 1, by = 2)]), 
                  sigdig)
    names(eps) <- as.character(x$nonmem$problem$estimation$epsshrink[seq(2, 
                                                                         length(x$nonmem$problem$estimation$epsshrink) - 1, by = 2)])
  } else {
    eta <- signif(as.numeric(x$nonmem$problem$estimation$etashrinksd[seq(1, 
                                                                         length(x$nonmem$problem$estimation$etashrinksd) - 1, by = 2)]), 
                  sigdig)
    names(eta) <- as.character(x$nonmem$problem$estimation$etashrinksd[seq(2, 
                                                                           length(x$nonmem$problem$estimation$etashrinksd) - 1, by = 2)])
    eps <- signif(as.numeric(x$nonmem$problem$estimation$epsshrinksd[seq(1, 
                                                                         length(x$nonmem$problem$estimation$epsshrinksd) - 1, by = 2)]), 
                  sigdig)
    names(eps) <- as.character(x$nonmem$problem$estimation$epsshrinksd[seq(2, 
                                                                           length(x$nonmem$problem$estimation$epsshrinksd) - 1, by = 2)]) 
    etavr <- signif(as.numeric(x$nonmem$problem$estimation$etashrinkvr[seq(1, 
                                                                           length(x$nonmem$problem$estimation$etashrinkvr) - 1, by = 2)]), 
                    sigdig)
    names(etavr) <- as.character(x$nonmem$problem$estimation$etashrinkvr[seq(2, 
                                                                             length(x$nonmem$problem$estimation$etashrinkvr) - 1, by = 2)])
    epsvr <- signif(as.numeric(x$nonmem$problem$estimation$epsshrinkvr[seq(1, 
                                                                           length(x$nonmem$problem$estimation$epsshrinkvr) - 1, by = 2)]), 
                    sigdig)
    names(epsvr) <- as.character(x$nonmem$problem$estimation$epsshrinkvr[seq(2, 
                                                                             length(x$nonmem$problem$estimation$epsshrinkvr) - 1, by = 2)]) 
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
    if(!is.null(x$nonmem$problem$estimation$etashrink)) {
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
