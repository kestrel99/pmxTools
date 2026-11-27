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
      stats::setNames(tabFiles, tabFiles)
    } else {
      stats::setNames(
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
    if (!quiet) {
      message(
        "The following files were not found and therefore cannot be loaded: ",
        paste0('"', tabFiles[!mask_found_files], '"', collapse=", ")
      )
    }
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
        bind_nmtables(loaded_files)
      } else if (length(lngths) == 2) {
        # There is likely a combination of firstonly and non-firstonly tables.
        # Per-record tables share NONMEM's record order; per-subject tables
        # need not share an ID order, so each is joined by ID separately.
        long <- bind_nmtables(loaded_files[file.df$filedim == max(lngths)])
        Reduce(f=merge_firstonly, x=loaded_files[file.df$filedim == min(lngths)], init=long)
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

# Column-bind tables with the same number of rows, which NONMEM writes in the
# same record order.  Columns repeated across tables (e.g. ID and TIME) are kept
# once, after checking that they agree row by row, i.e. that the rows line up.
bind_nmtables <- function(tables) {
  ret <- tables[[1]]
  for (tab in tables[-1]) {
    shared <- intersect(names(ret), names(tab))
    misaligned <- shared[!vapply(shared, function(col) same_values(ret[[col]], tab[[col]]), logical(1))]
    if (length(misaligned) > 0) {
      stop("Tables with the same number of rows disagree in shared column(s) ",
           paste(misaligned, collapse=", "), ", so their rows do not line up. ",
           "Use output_type = \"list\" to read them separately.", call.=FALSE)
    }
    ret <- cbind(ret, tab[, setdiff(names(tab), shared), drop=FALSE])
  }
  ret
}

# Element-wise equality allowing for the rounding of different $TABLE formats
same_values <- function(x, y, tolerance=1e-4) {
  both_na <- is.na(x) & is.na(y)
  close <- abs(x - y) <= pmax(tolerance * pmax(abs(x), abs(y)), 1e-12)
  all(both_na | (!is.na(close) & close))
}

# Merge a per-subject (FIRSTONLY) table onto per-record tables by ID.  Columns
# present in both are taken from the per-record tables, all per-record rows are
# kept in their original order, and subjects without a FIRSTONLY row get NA.
merge_firstonly <- function(long, short) {
  if (!("ID" %in% names(long) && "ID" %in% names(short))) {
    stop("Tables with different numbers of rows can only be merged by an `ID` column present in all of them. ",
         "Use output_type = \"list\" to read them separately.", call.=FALSE)
  }
  if (anyDuplicated(short$ID)) {
    stop("The shorter (FIRSTONLY) tables contain more than one row per `ID` and cannot be merged by subject. ",
         "Use output_type = \"list\" to read them separately.", call.=FALSE)
  }
  missing_ids <- setdiff(long$ID, short$ID)
  if (length(missing_ids) > 0) {
    warning("No FIRSTONLY table row for ID: ", paste(missing_ids, collapse=", "),
            "; their per-subject columns are NA.", call.=FALSE)
  }
  short <- short[, c("ID", setdiff(names(short), names(long))), drop=FALSE]
  row_order <- seq_len(nrow(long))
  ret <- merge(cbind(long, .row_order=row_order), short, by="ID", all.x=TRUE, sort=FALSE)
  ret <- ret[order(ret$.row_order), c(names(long), setdiff(names(short), "ID")), drop=FALSE]
  rownames(ret) <- NULL
  ret
}

# Keep only data rows (dropping repeated headers) and convert every column to
# numeric without changing the table's shape.
clean_nmtable <- function(tab) {
  automatic_rownames <- .row_names_info(tab) < 0L
  pattern <- "^[[:space:]]*([[:alpha:]].*)"
  tab <- tab[!grepl(pattern, tab[[1]]), , drop=FALSE]
  tab[] <- lapply(tab, as.numeric)
  # Renumber automatic row names after dropping rows; keep real labels (e.g.
  # parameter names read with row.names=1)
  if (automatic_rownames) rownames(tab) <- NULL
  tab
}

read_nmtable_single <- function(filename, quiet) {
  if (!quiet) {
    message("Reading ", filename)
  }
  
  ## Drop the "TABLE NO." banners. NONMEM writes one per table, so there are
  ## several (each followed by a repeated header row) when NSIM > 1.
  lines <- readLines(filename)
  banners <- grepl(pattern="^TABLE NO", x=lines)
  if (sum(banners) > 1 && !quiet) {
    message(
      filename, " contains ", sum(banners), " tables, as when running multiple ",
      "simulations in NONMEM (NSIM > 1). They will be stacked."
    )
  }
  lines <- lines[!banners]

  ## Check which type of separator we have in our tables
  sep.char <- if (grepl(pattern=",", x=lines[1])) "," else ""
  ret <- utils::read.table(text=lines, header=TRUE, sep=sep.char)
  # scan to ensure no remaining column headers or other text pollution
  clean_nmtable(ret)
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
    
    tmptab <- utils::read.table(
      textConnection(file_data[start_line:end_line]),
      header=header, ...
    )
    
    # scan to ensure no remaining column headers or other text pollution
    ret[[current_table_name]] <- clean_nmtable(tmptab)
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
