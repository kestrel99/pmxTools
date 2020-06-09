#' Read in the NONMEM variance-covariance matrix.
#'
#' @param fileName Root filename for the NONMEM run (e.g. "run315").
#'
#' This function reads the ".cov" NONMEM output table, and will return an error if this
#' is missing.
#'
#' @return A symmetrical variance-covariance matrix covering all model parameters.
#'
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' nmVcov <- read_nmcov("run315")
#' }
#' @import utils
#' @export

read_nmcov <- function(fileName) {
  if (!file.exists(fileName)) {
    fileName_orig <- fileName
    fileName <- paste0(fileName, ".cov")
    if (!file.exists(fileName)) {
      stop("Neither ", fileName, " nor ", fileName_orig, " exist.")
    } else {
      message(fileName_orig, " does not exist, trying ", fileName, ".")
    }
  }
  as.matrix(
    read.table(fileName, skip=1, header=TRUE, row.names=1)
  )
}
