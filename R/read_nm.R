#' Read NONMEM 7.2+ output into a list of lists.
#'
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
#'
#' @export
#' @importFrom xml2 read_xml as_list
read_nm <- function(fileName) {

  if(length(grep(".xml$", fileName))==0) {
    fileName <- paste(fileName, ".xml", sep="")
  }

  # nmFile <- XML::xmlTreeParse(fileName)
  # XML::xmlToList(nmFile)
  nmFile <- xml2::read_xml(fileName, ".xml", sep = "")
  xml2::as_list(nmFile)$output
}


