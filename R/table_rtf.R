#' Read NONMEM output into a list.
#'
#' \code{table_rtf} generates an RTF table from a data frame.
#'
#' @param df A data frame.
#' @param outFile A filename for writing the table to. If \code{NULL}, writes to console.
#' @param rtfFile If \code{TRUE} (the default), then add RTF tabs to generate a fully formatted RTF file.
#' @param boldHeader If \code{TRUE}, make the header bold.
#' @param rowNames If \code{TRUE}, include row names in the table. Default is \code{FALSE}.
#' @param ... Other formatting options for the table body.
#'
#' @return An RTF table based on the data frame provided.
#'
#' @examples
#' \dontrun{
#' scm <- read_scm("E:/DrugX/ModelDevelopment/scm310")
#' myRTF <- table_rtf(scm$forwardSummary)
#' }
#'
#' @references \url{http://www.r-bloggers.com/another-solution-to-the-r-to-word-table-problem/}
#'
#' @export

table_rtf <- function(df, outFile=NULL, rtfFile=TRUE, boldHeader=TRUE,
                     rowNames=FALSE, ...) {

  if (!is.null(outFile)) sink(outFile)

  df.nrow <- nrow(df)
  df.ncol <- ncol(df)

  if(rowNames==T) {
    rn <- 1
  } else {
    rn <- 0
  }

  if (rtfFile) {
    cat("{\\rtf1\n")
  }
  # populate header row
  cat("\\trowd\\trautofit1\\intbl\n")
  j <- 1
  for (i in 1:(df.ncol+rn)) {
    cat("\\cellx",j,'\n',sep='')
    j<-j+1
  }
  cat("{\n")
  # loop through and write column headers
  if (rowNames==T) cat(" \\cell\n")
  for (i in 1:df.ncol) {
    if (boldHeader) {
      cat('\\b ',colnames(df)[i],"\\b0\\cell \n",sep='')
    } else {
      cat(colnames(df)[i],"\\cell \n",sep='')
    }
  }
  cat("}\n")
  cat("{\n")
  cat("\\trowd\\trautofit1\\intbl\n")

  j<-1
  for (i in 1:(df.ncol+rn)) {
    cat("\\cellx",j,'\n',sep='')
    j<-j+1
  }
  cat("\\row }\n")

  # write table contents
  for (k in 1:df.nrow) {
    cat("\\trowd\\trautofit1\\intbl\n")

    j<-1
    for (i in 1:(df.ncol+rn)) {
      cat("\\cellx",j,'\n',sep='')
      j<-j+1
    }
    cat("{\n")
    if (rowNames==T) cat(rownames(df)[k],'\\cell\n',sep='')
    for (i in 1:df.ncol) {
      cat(format(df[k,i],...),"\\cell \n",sep='')
    }
    cat("}\n")
    cat("{\n")
    cat("\\trowd\\trautofit1\\intbl\n")
    j<-1
    for (i in 1:(df.ncol+rn)) {
      cat("\\cellx",j,'\n',sep='')
      j<-j+1
    }
    cat("\\row }\n")
  }
  if (rtfFile) {
    cat("}\n")
  }
  if (!is.null(outFile)) sink()
}

