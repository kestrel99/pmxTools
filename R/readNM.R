#' Read NONMEM 7.2+ output into a list of lists.
#'
#' @param fileName A NONMEM XML output file (e.g. "run315.xml").
#'
#' @return A list of lists corresponding to the NONMEM output object.
#'
#' @examples
#' nmOutput <- readNM("run315.xml")

readNM <- function(fileName) {

  if(length(grep(".xml$", fileName))==0) {
    fileName <- paste(fileName, ".xml", sep="")
  }

  nmFile <- xmlTreeParse(fileName)
  xmlToList(nmFile)

}
