#' Extract structural model parameter estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{readNM}}.
#'
#' @return A named vector of NONMEM model parameter estimates.
#'
#' @examples
#'  nmOutput <- readNM("run315.xml")
#'  thetas <- getTheta(nmOutput)

getTheta <- function(x) {

  out <- as.numeric(x$nonmem$problem$estimation$theta[1,])
  names(out) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")
  out
}
