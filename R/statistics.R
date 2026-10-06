#' Calculate geometric mean
#'
#' @param x Numeric vector.
#' @param na.rm Flag for removing \code{NA} values (defaults to \code{FALSE}).
#' @param neg.rm Flag for removing negative or zero values (defaults to \code{FALSE}).
#'  
#' @return The geometric mean. \code{NA} is returned if there are any non-positive elements in \code{x}.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#'
#' @examples
#' gm(c(0.5, 7, 8, 5))
#'
#' @export

gm <- function(x, na.rm=FALSE, neg.rm=FALSE) {
  z <- x
  if (na.rm) z <- z[!is.na(z)]
  if (neg.rm) z <- z[z > 0]

  if(any(is.na(z))) {
    out <- as.numeric(NA)
  }
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
#' @param na.rm A logical value indicating whether `NA` values should be stripped before the computation proceeds.
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
#' 
#' \url{http://onbiostatistics.blogspot.com/2008/07/geometric-statistics-geometric-cv-vs.html}
#' @export

gcv_convert <- function(gvar=gsd^2, gsd) {
  if (!xor(missing(gvar), missing(gsd))) {
    stop("Only one of `gvar` or `gsd` may be provided at a time.")
  }
  100*sqrt(exp(gvar)-1)
}

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


