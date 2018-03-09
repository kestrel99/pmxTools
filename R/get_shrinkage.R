#' Extract shrinkage estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param output A flag specifying the shrinkage estimates to be output. Valid flag
#' values are \code{eta} (the default), \code{epsilon}, or \code{all}.
#' @param sigdig Specifies the number of significant digits to be provided (default=3).
#'
#' @return A named vector of NONMEM shrinkage estimates, or in the case of \code{all},
#' a list of named vectors.
#'
#' \code{eta} returns a vector of ETA shrinkages, as reported by NONMEM.
#' \code{epsilon} returns EPSILON shrinkage, as reported by NONMEM.
#' \code{all} returns both ETA and EPSILON shrinkage estimates as a list of vectors.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  nmOutput <- read_nm("run315.xml")
#'  shr <- get_shrinkage(nmOutput, output="all")
#' }
#'
#' @export

get_shrinkage <- function(x, output = "eta", sigdig = 3) {
  if (!(output %in% c("eta", "epsilon", "all"))) {
    stop("Please select a valid output option (eta, epsilon, all).")
  }

  eta <-
    signif(as.numeric(x$nonmem$problem$estimation$etashrink[seq(1,
                                                                length(x$nonmem$problem$estimation$etashrink) -
                                                                  1, by = 2)]),
           sigdig)
  names(eta) <-
    as.character(x$nonmem$problem$estimation$etashrink[seq(2,
                                                           length(x$nonmem$problem$estimation$etashrink) -
                                                             1, by = 2)])

  eps <-
    signif(as.numeric(x$nonmem$problem$estimation$epsshrink[seq(1,
                                                                length(x$nonmem$problem$estimation$epsshrink) -
                                                                  1, by = 2)]),
           sigdig)
  names(eps) <-
    as.character(x$nonmem$problem$estimation$epsshrink[seq(2,
                                                           length(x$nonmem$problem$estimation$epsshrink) -
                                                             1, by = 2)])

  if (output == "eta")
    out <- eta
  if (output == "epsilon")
    out <- eps
  if (output == "all") {
    out <- list()
    out$eta <- eta
    out$epsilon <- eps
  }

  out
}
