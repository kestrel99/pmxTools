#' Calculate a geometric coefficient of variation.
#'
#' @param x A vector.
#' @param na.rm Flag for removing \code{NA} values (defaults to \code{FALSE}).
#' @param neg.rm Flag for removing negative or zero values (defaults to \code{FALSE}).
#'
#' @return The geometric coefficient of variation of the input vector. If \code{neg.rm} is \code{FALSE} and values <= 0 are present, \code{NA} will be returned.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  gcv(myvector)
#' }
#'
#' @export

gcv <- function(x, na.rm=F, neg.rm=F) {
  xx <- x
  if (na.rm) xx <- xx[!is.na(xx)]
  if (neg.rm) xx <- xx[xx > 0]
  
  if(any(xx<=0)) {
    out <- NA
  } else {
    out <- sqrt(exp(sd(log(xx[!is.na(xx)]))^2) - 1)
  }
  out
}

