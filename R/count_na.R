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
  length(x[is.na(x)])
}