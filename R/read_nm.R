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
  
  testFile <- readLines(fileName)
  if(testFile[length(testFile)] != "</nm:output>") {
    warning("Invalid XML: ", fileName)
    return(NULL)    
  } else {
    nmFile <- xml2::read_xml(fileName, ".xml", sep = "")
    xml2::as_list(nmFile)$output
  }
  
}
