#' Format a number with the correct number of significant digits and trailing zeroes.
#'
#' @param x A vector of numeric values.
#' @param digits The number of significant digits values should have (defaults to 3).
#' 
#' @return A string containing the properly-formatted number.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  fmt_signif(c(36.44, 0.0002, 3336.7), digits=3)
#' }
#'
#' @export

fmt_signif <- function(x, digits=3) {
  x[is.na(x)] <- -1000000000
  x <- signif(x, digits)
  om <- floor(log10(abs(x)))
  om[is.infinite(om)] <- 0
  dp <- digits - om - 1
  dp[which(dp<0)] <- 0
  
  x <- sprintf(paste("%.",dp, "f", sep=""), x)
  x[x=="-1000000000"] <- "-"
  return(x)
}