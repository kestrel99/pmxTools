#' Reads NONMEM output tables.
#'
#' @inheritParams read_nm_all
#' @param tableFiles NONMEM table files to be read.
#' @param tabSuffix Table file suffix.
#' @param tableNames List of root table names, using the Xpose naming convention
#'   as the default.
#' @param output_type Should output be a "data.frame" where all results are
#'   merged or a "list" of data.frames.
#'
#' @return A data.frame or list of data.frames depending on the
#'   \code{output_type} argument.
#'
#' @note Adapted from Xpose 4 (\url{https://CRAN.R-project.org/package=xpose4}).
#' @references NONMEM (\url{https://www.iconplc.com/solutions/technologies/nonmem})
#' @references Jonsson EN, Karlsson MO. Xpose--an S-PLUS based population
#'   pharmacokinetic/pharmacodynamic model building aid for NONMEM. Comput
#'   Methods Programs Biomed. 1999 Jan;58(1):51-64
#' @family NONMEM reading
#' @author Bill Denney, Justin Wilkins, Niclas Jonsson, Andrew Hooker
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
      stop("Run number must be specified if no table files are provided.")
      return(NULL)
    }
    tabFiles <- paste0(tableNames, runNo, tabSuffix)
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
  # scan to ensure no remaining column headers or other text pollution

  # Define the regular expression pattern
  pattern <- "^[[:space:]]*([[:alpha:]].*)"
  
  # Find the indices of rows that match the pattern
  indices <- which(!grepl(pattern, ret[[1]]))
  ret <- ret[indices,]
  ret <- as.data.frame(sapply(ret, as.numeric))
  
  ret
}

#' Read (single or) multiple NONMEM tables from a single file
#' 
#' @param fileName The filename to read from
#' @param header,... Arguments passed to read.table
#' @param simplify If a single table is present, return a data.frame instead of
#'   a list of data.frames?
#' @param table_start_pattern What should be found to start a new table?
#' @return A list of data.frames, or if only one is present and simplify=TRUE, a
#'   data.frame.
#' @family NONMEM reading
#' @author Bill Denney
#' @examples
#' \dontrun{
#' read_nm_multi_table("run1.cov", row.names=1)
#' }
#' @export
read_nm_multi_table <- function(fileName, header=TRUE, ..., simplify=TRUE, table_start_pattern="^TABLE NO") {
  file_data <- readLines(fileName)
  new_tables <- grep(x=file_data, pattern=table_start_pattern)
  ret <- list()
  for (idx in seq_along(new_tables)) {
    start_line <- new_tables[idx] + 1
    end_line <- 
      if (idx == length(new_tables)) {
        length(file_data)
      } else {
        new_tables[idx + 1] - 1
      }
    current_table_name <- file_data[new_tables[idx]]
    
    tmptab <- read.table(
      textConnection(file_data[start_line:end_line]),
      header=header, ...
    )
    
    # scan to ensure no remaining column headers or other text pollution
    
    # Define the regular expression pattern
    pattern <- "^[[:space:]]*([[:alpha:]].*)"
    
    # Find the indices of rows that match the pattern
    indices <- which(!grepl(pattern, tmptab[[1]]))
    tmptab <- tmptab[indices,]
    tmptab <- as.data.frame(sapply(tmptab, as.numeric))    
    ret[[current_table_name]] <- tmptab
  }
  if (simplify & (length(ret) == 1)) {
    # Simplify to just the matrix in the common case of a single estimation
    # step
    ret <- ret[[1]]
  }
  ret
}

#' Read a standard NONMEM extension file
#' 
#' @param fileName The filename (with directory name, if applicable) to read
#'   (with or without the extension)
#' @param extension The file extension to optionally append (preferably starting
#'   with a ".")
#' @param directory The directory to look for files within.  If NULL, uses the
#'   current directory.
#' @param ... Passed to \code{read_nm_multi_table()}
#' @return NULL if the file does not exist or the value of
#'   \code{read_nm_multi_table()} if it does exist.
#' @examples
#' \dontrun{
#' read_nm_std_ext("run1", "phi")
#' }
#' @export
read_nm_std_ext <- function(fileName, extension, directory=NULL, ...) {
  if (!is.character(fileName)) fileName <- as.character(fileName)
  if (!startsWith(extension, prefix=".")) {
    extension <- paste0(".", extension)
  }
  fileName_read <- check_file_exists(fileName=fileName, ext=extension, directory=directory)
  if (is.null(fileName_read)) {
    warning("Could not find file: ", fileName, ", with extension: ", extension)
    return(NULL)
  }
  read_nm_multi_table(fileName=fileName_read, ...)
}
