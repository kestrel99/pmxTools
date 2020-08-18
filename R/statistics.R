#' Calculate geometric mean
#'
#' @param x Numeric vector.
#'  
#' @return The geometric mean. \code{NA} is returned if there are any non-positive elements in \code{x}.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#'
#' @examples
#' gm(c(0.5, 7, 8, 5))
#'
#' @export

gm <- function(x) {
  z <- x
  if(length(z[z<=0])>0) {
    out <- NA
  } else {
    out <- exp(mean(log(z)))
  }
  out
}

#' Calculate percentage coefficient of variation
#'
#' @param x Numeric vector.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#' 
#' @return The percentage coefficient of variation.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#'
#' @examples
#' pcv(rnorm(50, 5, 7.56))
#'
#' @export

pcv <- function(x, na.rm=FALSE) {
  sd(x, na.rm=na.rm)/mean(x, na.rm=na.rm)*100 
}

#' Convert geometric variance or standard deviation to a geometric coefficient of variation
#'
#' The equation used is: \code{100*sqrt(exp(gvar)-1)}
#'
#' @param gvar The geometric variance (note that this is the variance not a
#'   vector of values to compute the gcv from)
#' @param gsd The geometric standard deviation
#'
#' @return Geometric coefficient of variation
#' 
#' @author Bill Denney
#'
#' @examples
#' gcv_convert(0.2)
#' gcv_convert(gsd=0.2)
#' @references 
#' \url{https://www.cognigen.com/nonmem/nm/99feb042003.html}
#' 
#' \url{http://onbiostatistics.blogspot.com/2008/07/geometric-statistics-geometric-cv-vs.html}
#' @export
gcv_convert <- function(gvar=gsd^2, gsd) {
  if (!xor(missing(gvar), missing(gsd))) {
    stop("Only one of `gvar` or `gsd` may be provided at a time.")
  }
  100*sqrt(exp(gvar)-1)
}
