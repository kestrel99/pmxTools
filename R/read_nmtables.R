#' Reads NONMEM output tables.
#'
#' @param tableFiles NONMEM table files to be read.
#' @param runNo Run number.
#' @param tabSuffix Table file suffix.
#' @param tableNames List of root table names, using the Xpose naming convention
#'   as the default.
#' @param quiet Flag for displaying intermediate output (defaults to
#'   \code{FALSE}).
#' @param directory The directory to look for files within.  If NULL, uses the
#'   current directory.
#' @param ... Additional arguments.
#'
#' @return A data frame.
#'
#' @note Adapted from Xpose 4 (\url{https://CRAN.R-project.org/package=xpose4}).
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' @seealso Jonsson EN, Karlsson MO. Xpose--an S-PLUS based population
#'   pharmacokinetic/pharmacodynamic model building aid for NONMEM. Comput
#'   Methods Programs Biomed. 1999 Jan;58(1):51-64
#' @author Bill Denney, Justin Wilkins, Niclas Jonsson, Andrew Hooker
#'
#' @examples
#' \dontrun{
#' tables <- read_nmtables(runNo=315)
#' }
#'
#' @import utils
#' @export

read_nmtables <- function(tableFiles = NULL,
                          runNo      = NULL,
                          tabSuffix  = "",
                          tableNames = c("sdtab", "mutab", "patab", "catab",
                                         "cotab", "mytab", "extra", "xptab"),
                          quiet = FALSE,
                          directory = NULL,
                          output_type = c("data.frame", "list"),
                          ...) {
  output_type <- match.arg(output_type)

  # Determine the table file name
  if (is.null(tableFiles)){
    if(is.null(runNo)) {
      stop("Run number must be specified if no table files are provided.\n")
      return(NULL)
    }
    tabFiles <- sapply(tableNames, paste, runNo, tabSuffix, sep="")
  } else {
    tabFiles <- tableFiles
  }
    
  tabFiles <-
    if (is.null(directory)) {
      setNames(tabFiles, tabFiles)
    } else {
      setNames(
        file.path(directory, tabFiles),
        tabFiles
      )
    }
  mask_found_files <- sapply(X=tabFiles, FUN=file.exists)
  if (all(!mask_found_files)) {
    warning(
      "No files to load were found.  The following files do not exist: ",
      paste0('"', tabFiles[!mask_found_files], '"', collapse=", ")
    )
    return(NULL)
  } else if (any(!mask_found_files)) {
    message(
      "The following files were not found and therefore cannot be loaded: ",
      paste0('"', tabFiles[!mask_found_files], '"', collapse=", ")
    )
    tabFiles <- tabFiles[mask_found_files]
  }

  ## Read in the table files
  loaded_files <- 
    lapply(
      X=tabFiles,
      FUN=read_nmtable_single,
      quiet=quiet
    )

  if (output_type %in% "data.frame") {
    ## Check if the files have the same length
    file.df <-
      data.frame(
        seen.files=names(loaded_files),
        filedim=sapply(X=loaded_files, FUN=nrow),
        stringsAsFactors=FALSE
      )
    lngths  <- unique(file.df$filedim)
    
    ret <-
      if (length(lngths) == 1) {
        Reduce(f=cbind, x=loaded_files)
      } else if (length(lngths) == 2) {
        # There is likely a combination of firstonly and non-firstonly tables.
        short <- Reduce(f=cbind, x=loaded_files[lngths == min(lngths)])
        long <- Reduce(f=cbind, x=loaded_files[lngths == max(lngths)])
        merge(long, short)
      } else {
        message(
          "The table files associated with this run number (", runNo, ") appear ",
          "to have more than two different lengths. ",
          "Please check your output, it is likely files have been modified or the ",
          "$TABLE step has failed or request 'list' output_type."
        )
        NULL
      }
  } else {
    ret <- loaded_files
  }
  ret
}

read_nmtable_single <- function(filename, quiet) {
  if (!quiet) {
    message("Reading ", filename)
  }
  
  ## Check which type of separator we have in our tables
  header.line <- scan(file=filename, nlines=1, skip=1, what="character", sep="\n", quiet=TRUE)
  sep.char <- ""
  if (length(grep(pattern=",", x=header.line)) != 0) {
    sep.char <- ","
  }

  ## Check if we have unequal number of fields in the file
  ## used for multiple simulations
  fields.per.line      <- count.fields(filename)
  fields.in.first.line <- fields.per.line[1]
  fields.in.rest       <- fields.per.line[-1]
  if((length(unique(fields.in.rest)) != 1) ||
     (all(fields.in.first.line == fields.in.rest))) {
    if (!quiet) {
      message(
        filename, " contains varying numbers of fields.  ",
        "This may be due to multiple TABLE and header rows ",
        "caused by running multiple simulations in NONMEM (NSIM > 1).",
        "Will attempt to remove these rows. Please be patient..."
      )
    }
    tmp   <- readLines(filename, n=-1)
    inds  <- grep(pattern="TABLE", x=tmp)
    if (length(inds) != 1) {
      inds  <- inds[c(2:length(inds))]
      inds2 <- inds+1
      ret <-
        read.table(
          textConnection(tmp[-c(inds, inds2)]),
          skip=2, header=TRUE, sep=sep.char
        )
    } else {
      ret <-
        read.table(filename, skip=2, header=TRUE, sep=sep.char)
    }
  } else {
    ret <-
      read.table(filename, skip=1, header=TRUE, sep=sep.char)
  }
  ret
}
