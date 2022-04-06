#' Read NONMEM output into a list.
#'
#' \code{read_nmext} returns a summary of a given NONMEM run, including
#' termination messages, parameter estimates, and precision estimates.
#' Minimally, the NONMEM output and '.ext' files must be available.
#'
#' @inheritParams read_nm_all
#' @param fileName A NONMEM output file prefix, without extension (e.g.
#'   "run315").
#' @param fileExt  The file extension for NONMEM output, set to ".lst" by
#'   default.
#' @param estNo The estimation number to report (by default, if only one
#'   estimation step is present, that will be reported; if multiple are
#'   reported, the last will be reported by default).
#'
#' @return A list of lists, containing 'Termination' (summary of NONMEM's
#'   termination output, including shrinkages and ETABAR estimates), 'OFV' (the
#'   objective function value), 'Thetas' (a vector of structural parameter
#'   estimates, or THETAs), 'Omega', a list of lists containing the OMEGA
#'   matrix, 'Sigma', a list of lists containing the SIGMA matrix, 'seThetas', a
#'   vector of standard errors for THETAs, 'seOmega', a list of lists containing
#'   standard errors for the OMEGA matrix, and 'seSigma', a list of lists
#'   containing standard errors for the SIGMA matrix.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' @family NONMEM reading
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' read_nmext("run315")
#' read_nmext("run315", ".nmlst")
#' }
#' @export
read_nmext <- function(fileName, fileExt = ".lst", directory=NULL, quiet=FALSE, estNo=NULL, ...) {
  fileName_read <- check_file_exists(fileName, fileExt, directory=directory)
  if (is.null(fileName_read)) {
    warning("Could not find file: ", fileName)
    return(NULL)
  }
  if (!quiet) {
    message("Reading ", fileName_read)
  }
  nmFile <-
    scan(
      fileName_read,
      sep = "\n",
      what = character(),
      quiet = TRUE
    )
  ret <- parse_nmext(nmFile, fileName_read, estNo=estNo)
  ret$raw_lst <- nmFile
  ret
}

parse_nmext <- function(nmFile, fileName_read, estNo=NULL) {
  minStart <- grep(pattern="#TERM:", x=nmFile)
  minEnd   <- grep(pattern="#TERE:", x=nmFile)
  if (!is.null(estNo)) {
    if (estNo > length(minStart)) {
      stop("Requested estimation number, ", estNo, " does not exist or did not complete successfully.")
    }
    minStart <- minStart[estNo]
    minEnd <- minEnd[estNo]
  } else if (length(minStart) > 1) {
    warning("More than one set of termination messages in this file. Returning information for the final estimation.")
    minStart <- max(minStart)
    minEnd <- max(minEnd)
  }
  if (length(minStart) == 0 || length(minEnd) == 0) {
    termMsg <- NULL
  } else {
    termMsg <- nmFile[(minStart + 1):(minEnd - 1)]
    termMsg <- substring(termMsg, 2)
  }
  extFileName <- paste0(sub("\\.\\w*$", "", fileName_read), ".ext")
  if (!file.exists(extFileName)) {
    warning(paste(
      "Could not find the raw results file (",
      extFileName,
      ") for  ",
      fileName_read,
      "."
    ))
    return(list(NULL))
  } else {
    extData <- read_nm_multi_table(extFileName, header = TRUE, simplify=FALSE)
    # Regardless, we need to get it out of a list
    extData <-
      if (is.null(estNo)) {
        extData[[length(extData)]]
      } else if (estNo <= length(extData)) {
        extData[[estNo]]
      } else {
        stop("Fewer tables than the requested estimation number (", estNo, ") in ", extFileName)
      }
  }

  ofv          <- extData$OBJ[extData$ITERATION == -1e+09]
  finalEstLine <-
    as.numeric(extData[extData$ITERATION == -1e+09,])
  thetaList    <- finalEstLine[grep("THETA", names(extData))]
  omegaCount   <- grep("OMEGA", names(extData))
  omegaList    <- list()
  seenOmega    <- NULL
  seenOmega[1:100] <- 0

  for (omegaNo in omegaCount) {
    omegaName <- names(extData)[omegaNo]
    omegaCol  <-
      as.numeric(sub("OMEGA\\.\\w*\\.", "", omegaName, perl = TRUE))
    omegaRow  <-
      as.numeric(sub(
        "\\.\\w*\\.$",
        "",
        sub("OMEGA\\.", "", omegaName, perl = TRUE),
        perl = TRUE
      ))
    seenOmega[omegaRow] <- seenOmega[omegaRow] + 1
    if (seenOmega[omegaRow] == 1) {
      omegaList[omegaRow][omegaCol] <- finalEstLine[omegaNo]
    }
    else {
      omegaList[[omegaRow]] <-
        c(omegaList[[omegaRow]], as.numeric(finalEstLine[omegaNo]))
    }
  }

  sigmaCount <- grep("SIGMA", names(extData))
  sigmaList  <- list()
  seenSigma  <- NULL
  seenSigma[1:100] <- 0
  if (length(sigmaCount) != 0) {
    for (sigmaNo in sigmaCount) {
      sigmaName <- names(extData)[sigmaNo]
      sigmaCol <-
        as.numeric(sub("SIGMA\\.\\w*\\.", "", sigmaName, perl = TRUE))
      sigmaRow <-
        as.numeric(sub(
          "\\.\\w*\\.$",
          "",
          sub("SIGMA\\.", "", sigmaName, perl = TRUE),
          perl = TRUE
        ))
      seenSigma[sigmaRow] <- seenSigma[sigmaRow] + 1
      if (seenSigma[sigmaRow] == 1) {
        sigmaList[sigmaRow][sigmaCol] <- finalEstLine[sigmaNo]
      }
      else {
        sigmaList[[sigmaRow]] <-
          c(sigmaList[[sigmaRow]], as.numeric(finalEstLine[sigmaNo]))
      }
    }
  } else {
    sigmaList <- NULL
  }

  seEstLine <-
    as.numeric(extData[extData$ITERATION == -1000000001,])
  if (length(seEstLine)[1] == 0) {
    seThetas <- NULL
    seOmegas <- NULL
    seSigmas <- NULL
  } else {
    seThetas <- seEstLine[grep("THETA", names(extData))]
    seThetas[seThetas == 1e+10] <- NA

    omegaCount  <- grep("OMEGA", names(extData))
    seOmegaList <- list()
    seenSeOmega <- NULL
    seenSeOmega[1:100] <- 0
    for (omegaNo in omegaCount) {
      omegaName <- names(extData)[omegaNo]
      omegaCol <-
        as.numeric(sub("OMEGA\\.\\w*\\.", "", omegaName, perl = TRUE))
      omegaRow <-
        as.numeric(sub(
          "\\.\\w*\\.$",
          "",
          sub("OMEGA\\.", "", omegaName, perl = TRUE),
          perl = TRUE
        ))
      if (!is.na(seEstLine[omegaNo])) {
        if (seEstLine[omegaNo] == 1e+10)
          seEstLine[omegaNo] <- NA
      }
      seenSeOmega[omegaRow] <- seenSeOmega[omegaRow] + 1
      if (seenSeOmega[omegaRow] == 1) {
        seOmegaList[omegaRow][omegaCol] <- seEstLine[omegaNo]
      } else {
        seOmegaList[[omegaRow]] <-
          c(seOmegaList[[omegaRow]], as.numeric(seEstLine[omegaNo]))
      }
    }

    sigmaCount  <- grep("SIGMA", names(extData))
    seSigmaList <- list()
    seenSeSigma <- NULL
    seenSeSigma[1:100] <- 0
    if (length(sigmaCount) != 0) {
      for (sigmaNo in sigmaCount) {
        sigmaName <- names(extData)[sigmaNo]
        sigmaCol  <-
          as.numeric(sub("SIGMA\\.\\w*\\.", "", sigmaName, perl = TRUE))
        sigmaRow  <-
          as.numeric(sub(
            "\\.\\w*\\.$",
            "",
            sub("SIGMA\\.", "", sigmaName, perl = TRUE),
            perl = TRUE
          ))
        if (!is.na(seEstLine[sigmaNo])) {
          if (seEstLine[sigmaNo] == 1e+10)
            seEstLine[sigmaNo] <- NA
        }
        seenSeSigma[sigmaRow] <- seenSeSigma[sigmaRow] + 1
        if (seenSeSigma[sigmaRow] == 1) {
          seSigmaList[sigmaRow][sigmaCol] <- seEstLine[sigmaNo]
        } else {
          seSigmaList[[sigmaRow]] <-
            c(seSigmaList[[sigmaRow]], as.numeric(seEstLine[sigmaNo]))
        }
      }
    } else {
      seSigmaList <- NULL
    }
  }

  out <- list(
    Termination = termMsg,
    OFV = ofv,
    Thetas = thetaList,
    Omega = omegaList,
    Sigma = sigmaList,
    seThetas = seThetas,
    seOmega = seOmegaList,
    seSigma = seSigmaList
  )
  return(out)
}
