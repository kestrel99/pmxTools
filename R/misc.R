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