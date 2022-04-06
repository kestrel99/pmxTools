#' Calculate a geometric mean.
#'
#' @param x A vector.
#' @param na.rm Flag for removing \code{NA} values (defaults to \code{FALSE}).
#' @param neg.rm Flag for removing negative or zero values (defaults to \code{FALSE}).
#'
#' @return The geometric mean of the input vector. If \code{neg.rm} is \code{FALSE} and values <= 0 are present, \code{NA} will be returned.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  geomean(myvector)
#' }
#'
#' @export
#' 
geomean <- function(x, na.rm=F, neg.rm=F) {
  xx <- x
  if (na.rm) xx <- xx[!is.na(xx)]
  if (neg.rm) xx <- xx[xx > 0]
  
  if(any(xx<=0)) {
    out <- NA 
  } else {
    out <- exp(mean(log(xx))) 
  }
  out
}