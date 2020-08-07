#' Read all NONMEM files for a single NONMEM run.
#' 
#' @details The filename for loading is constructed as \code{paste(run_prefix,
#'   runNo)}.  To load a nonstandard file, simply set one of those values to
#'   \code{NULL}.
#' 
#' @param run_prefix The start to the name of the run.
#' @param runNo Run number.
#' @param directory The directory to look for files within.  If NULL, uses the
#'   current directory.
#' @param quiet Flag for displaying intermediate output.
#' @param ... Passed to each of the read functions (ignored in the functions).
#' @family NONMEM reading
#' @export
read_nm_all <- function(runNo, run_prefix="run", directory=NULL, quiet=FALSE, ...) {
  list(
    xml=read_nm(fileName=paste0(run_prefix, runNo), quiet=quiet, directory=directory, ...),
    ext=read_nmext(fileName=paste0(run_prefix, runNo), quiet=quiet, directory=directory, ...),
    cov=read_nmcov(fileName=paste0(run_prefix, runNo), quiet=quiet, directory=directory, ...),
    extra_files=
      lapply(
        X=setNames(nm=c("ext", "phi", "ets", "phm")),
        FUN=read_nm_std_ext,
        fileName=paste0(run_prefix, runNo),
        directory=directory
      ),
    tables=read_nmtables(runNo=runNo, directory=directory, ...)
  )
}

#' Check if a filename exists with or without the extension added.
#' 
#' @param fileName The filename that may exist
#' @param ext The filename (with dot, if applicable)
#' @return The filename if it exists, the filename with extension if filename
#'   does not exist and the filename with extension does exist, and NULL if
#'   neither exists.
#' @examples
#' check_file_exits("foo", ".bar")
#' @noRd
check_file_exists <- function(fileName, ext=NULL, directory=NULL) {
  if (!is.null(directory)) {
    fileName <- file.path(directory, fileName)
  }
  file_with_ext <- paste0(fileName, ext)
  if (file.exists(fileName)) {
    fileName
  } else if (file.exists(file_with_ext)) {
    file_with_ext
  } else {
    NULL
  }
}
