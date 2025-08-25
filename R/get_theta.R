#' Extract structural model parameter estimates and associated information
#' from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param output A flag specifying the matrix or matrices to be output. Valid flag values are \code{est} (the default),
#'  \code{se}, \code{rse}, \code{95ci}, or \code{all}.
#' @param sigdig Specifies the number of significant digits to be provided (default=6).
#' @param sep Specifies the separator character to use for 95% confidence intervals (default="-").
#' @param est.step Specifies which estimation step to return parameters from (default is the last).
#'
#' @return A named vector of NONMEM model parameter estimates, or in the case of \code{all},
#' a list of named vectors.
#'
#' @details Output options are as follows:
#' _est_ returns a vector of `THETA` values.
#' _se_ returns a vector of `THETA` standard errors.
#' _rse_ returns a vector of `THETA` relative standard errors (`se/est*100`).
#' _95ci_ returns a vector of the asymptotic 95% confidence intervals for the elements of `THETA` (`est +/- 1.96*se`).
#' _all_ returns all available `THETA` information as a list of named vectors.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  nmOutput <- read_nm("run315.xml")
#'  thetas <- get_theta(nmOutput)
#' }
#'
#' @export

get_theta <- function(x, output="est", sigdig=6, sep="-", est.step=NULL) {

  if (!(output %in% c("est", "se", "rse", "95ci", "all"))) {
    stop("Please select a valid output option (est, se, rse, 95ci, all).")
  }

  if(is.null(est.step)) {
    no_steps <- sum(stringr::str_count(names(x$nonmem$problem), "estimation"))
  } else {
    no_steps <- est.step
  }
  
  ind_est  <- match("estimation", names(x$nonmem$problem))-1+no_steps
  
  if(output=="est") {
    out <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$theta)), sigdig)
    names(out) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")
  }

  if(output=="se") {
    out <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$thetase)), sigdig)
    names(out) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")
  }

  if(output=="rse") {
    est <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$theta))
    se  <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$thetase))
    out <- signif(abs(se/est*100), sigdig)
    names(out) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")
  }

  if(output=="95ci") {
    est <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$theta))
    se  <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$thetase))
    outup <- fmt_signif(est + 1.96*se, digits=sigdig)
    outlo <- fmt_signif(est - 1.96*se, digits=sigdig)
    out <- paste(outlo, outup, sep=sep)
    names(out) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")
  }

  if(output=="all") {
    out <- list()

    out$Theta <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$theta)), sigdig)
    names(out$Theta) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")

    out$ThetaSE <- signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$thetase)), sigdig)
    names(out$ThetaSE) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")

    est <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$theta))
    se  <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$thetase))
    out$ThetaRSE <- signif(abs(se/est*100), sigdig)
    names(out$ThetaRSE) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")

    est <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$theta))
    se  <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$thetase))
    outup <- fmt_signif(est + 1.96*se, digits=sigdig)
    outlo <- fmt_signif(est - 1.96*se, digits=sigdig)
    out$Theta95CI <- paste(outlo, outup, sep=sep)
    names(out$Theta95CI) <- paste("THETA", 1:length(x$nonmem$problem[[ind_est]]$theta), sep="")
  }

  out
}
