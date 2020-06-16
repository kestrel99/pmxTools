#' Read NONMEM 7.2+ output into a list of lists.
#'
#' @inheritParams read_nm_all
#' @param fileName A NONMEM XML output file (e.g. "run315.xml").
#'
#' @return A list of lists corresponding to a NONMEM output object.
#' 
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' nmOutput <- read_nm("run315.xml")
#' }
#' @family NONMEM reading
#' @export
#' @importFrom XML xmlToList xmlTreeParse
read_nm <- function(fileName, directory=NULL, quiet=FALSE, ...) {
  fileName_read <- check_file_exists(filename=fileName, ext=".xml", directory=directory)
  if (is.null(fileName_read)) {
    warning("Could not find file: ", fileName)
    return(NULL)
  }
  if (!quiet) {
    message("Reading ", fileName_read)
  }
  nmFile <- XML::xmlTreeParse(fileName_read)
  XML::xmlToList(nmFile)
}
