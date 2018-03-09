#' Read in the NONMEM variance-covariance matrix.
#'
#' @param fileName Root filename for the NONMEM run (e.g. "run315").
#'
#' This function reads the ".cov" NONMEM output table, and will return an error if this
#' is missing.
#'
#' @return A symmetrical variance-covariance matrix covering all model parameters.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' nmVcov <- read_nmcov("run315")
#' }
#' @import utils
#' @export

read_nmcov <- function(fileName) {
  if(file.exists(paste(fileName, ".cov", sep=""))) {
    as.matrix(read.table(paste(fileName, ".cov", sep=""), skip=1, header=T, row.names=1))
  } else {
    stop(paste("File ", paste(fileName, ".cov", sep=""), " not found.", sep=""))
  }
}

