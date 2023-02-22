#' Count the number of NA values in a vector.
#'
#' @param x A vector.
#'
#' @return An integer containing the number of NA values in the input vector.
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  count_na(c(0,5,7,NA,3,3,NA))
#' }
#'
#' @export

count_na <- function(x) {
  n_na  <- length(x[is.na(x)])
  n_nan <- length(x[is.nan(x)])
  if(n_nan==1) warning(paste(n_nan, "NaN value included in the NA count."))
  if(n_nan>1) warning(paste(n_nan, "NaN values included in the NA count."))
  n_na
}