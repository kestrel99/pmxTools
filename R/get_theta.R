#' Extract structural model parameter estimates and associated information
#' from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param output A flag specifying the matrix or matrices to be output. Valid flag values are \code{est} (the default),
#'  \code{se}, \code{rse}, \code{95ci}, or \code{all}.
#' @param sigdig Specifies the number of significant digits to be provided (default=6).
#' @param sep Specifies the separator character to use for 95\% confidence intervals (default="-").
#'
#' @return A named vector of NONMEM model parameter estimates, or in the case of \code{all},
#' a list of named vectors.
#'
#' \code{est} returns a vector of THETA values.
#' \code{se} returns a vector of THETA standard errors.
#' \code{rse} returns a vector of THETA relative standard errors (se/est*100).
#' \code{95ci} returns a vector of the asymptotic 95\% confidence intervals for the elements of THETA (est +/- 1.96*se).
#' \code{all} returns all available THETA information as a list of named vectors.
#'
#' @examples
#' \dontrun{
#'  nmOutput <- read_nm("run315.xml")
#'  thetas <- get_theta(nmOutput)
#' }
#'
#' @export

get_theta <- function(x, output="est", sigdig=6, sep="-") {

  if (!(output %in% c("est", "se", "rse", "95ci", "all"))) {
    stop("Please select a valid output option (est, se, rse, 95ci, all).")
  }

  if(output=="est") {
    out <- signif(as.numeric(x$nonmem$problem$estimation$theta[1,]), sigdig)
    names(out) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")
  }

  if(output=="se") {
    out <- signif(as.numeric(x$nonmem$problem$estimation$thetase[1,]), sigdig)
    names(out) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")
  }

  if(output=="rse") {
    est <- as.numeric(x$nonmem$problem$estimation$theta[1,])
    se  <- as.numeric(x$nonmem$problem$estimation$thetase[1,])
    out <- signif(abs(se/est*100), sigdig)
    names(out) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")
  }

  if(output=="95ci") {
    est <- as.numeric(x$nonmem$problem$estimation$theta[1,])
    se  <- as.numeric(x$nonmem$problem$estimation$thetase[1,])
    outup <- signif(est + 1.96*se, sigdig)
    outlo <- signif(est - 1.96*se, sigdig)
    out <- paste(outlo, outup, sep=sep)
    names(out) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")
  }

  if(output=="all") {
    out <- list()

    out$Theta <- signif(as.numeric(x$nonmem$problem$estimation$theta[1,]), sigdig)
    names(out$Theta) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")

    out$ThetaSE <- signif(as.numeric(x$nonmem$problem$estimation$thetase[1,]), sigdig)
    names(out$ThetaSE) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")

    est <- as.numeric(x$nonmem$problem$estimation$theta[1,])
    se  <- as.numeric(x$nonmem$problem$estimation$thetase[1,])
    out$ThetaRSE <- signif(abs(se/est*100), sigdig)
    names(out$ThetaRSE) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")

    est <- as.numeric(x$nonmem$problem$estimation$theta[1,])
    se  <- as.numeric(x$nonmem$problem$estimation$thetase[1,])
    outup <- signif(est + 1.96*se, sigdig)
    outlo <- signif(est - 1.96*se, sigdig)
    out$Theta95CI <- paste(outlo, outup, sep=sep)
    names(out$Theta95CI) <- paste("THETA", x$nonmem$problem$estimation$theta[2,], sep="")
  }

  out
}
