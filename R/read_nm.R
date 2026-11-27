#' Read NONMEM 7.2+ output into a list of lists.
#'
#' @inheritParams read_nm_all
#' @param fileName A NONMEM XML output file (e.g. "run315.xml").
#'
#' @return A list of lists corresponding to a NONMEM output object.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/solutions/technologies/nonmem})
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' nmOutput <- read_nm("run315.xml")
#' }
#' @family NONMEM reading
#' @export
#' @importFrom xml2 read_xml as_list
read_nm <- function(fileName, directory=NULL, quiet=FALSE, ...) {

  fileName_read <- check_file_exists(fileName=fileName, ext=".xml", directory=directory)

  if (is.null(fileName_read)) {
    warning("Could not find file: ", fileName)
    return(NULL)
  }
  
  if (!quiet) {
    message("Reading ", fileName_read)
  }
  
  # Parse the document itself rather than checking its last line: a truncated
  # file (e.g. from a run still in progress) fails to parse, while trailing
  # whitespace and either line ending are valid XML.
  nmFile <- tryCatch(
    xml2::read_xml(fileName_read),
    error = function(e) {
      warning("Invalid XML: ", fileName_read, " (", conditionMessage(e), ")", call. = FALSE)
      NULL
    }
  )
  if (is.null(nmFile)) {
    return(NULL)
  }
  if (xml2::xml_name(nmFile) != "output") {
    warning("Invalid XML: ", fileName_read, " (root element is not <nm:output>)", call. = FALSE)
    return(NULL)
  }
  xml2::as_list(nmFile)$output
}
