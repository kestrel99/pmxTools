#' Read NONMEM 7.2+ output into a list of lists.
#'
#' @param fileName A NONMEM XML output file (e.g. "run315.xml").
#'
#' @return A list of lists corresponding to the NONMEM output object.
#'
#' @examples
#' \dontrun{
#' nmOutput <- read_nm("run315.xml")
#' }
#'
#' @export

read_nm <- function(fileName) {

  if(length(grep(".xml$", fileName))==0) {
    fileName <- paste(fileName, ".xml", sep="")
  }

  nmFile <- XML::xmlTreeParse(fileName)
  XML::xmlToList(nmFile)

}


