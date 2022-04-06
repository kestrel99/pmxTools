#' Read in the NONMEM variance-covariance matrix.
#'
#' @inheritParams read_nm_all
#' @param fileName Root filename for the NONMEM run (e.g. "run315").
#'
#' This function reads the ".cov" NONMEM output table, and will return an error if this
#' is missing.
#'
#' @return A symmetrical variance-covariance matrix covering all model parameters.
#'
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' nmVcov <- read_nmcov("run315")
#' }
#' @import utils
#' @family NONMEM reading
#' @export
read_nmcov <- function(fileName, quiet=FALSE, directory=NULL, ...) {
  fileName_read <- check_file_exists(fileName, ".cov", directory=directory)
  if (is.null(fileName_read)) {
    warning("Could not find file: ", fileName)
    return(NULL)
  } else {
    if (!quiet) {
      message("Reading ", fileName_read)
    }
    ret <-
      lapply(
        X=read_nm_multi_table(fileName=fileName_read, row.names=1, simplify=FALSE),
        FUN=as.matrix
      )
    if (length(ret) == 1) {
      ret <- ret[[1]]
    }
    ret
  }
}
