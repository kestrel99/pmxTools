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
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' @family NONMEM reading
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' read_nmext("run315")
#' read_nmext("run315", ".nmlst")
#' }
#' @export
read_nmext <- function(fileName, fileExt = ".lst", directory=NULL, quiet=FALSE, ...) {
  fileName_read <- check_file_exists(fileName, fileExt, directory=directory)
  if (is.null(fileName_read)) {
    warning("Could not find file: ", fileName)
    return(NULL)
  }
  if (!quiet) {
    message("Reading ", fileName_read)
  }
  nmFile <-
    scan(fileName_read,
         sep = "\n",
         what = character(),
         quiet = TRUE)
  minStart <- grep("#TERM:", nmFile)
  minEnd   <- grep("#TERE:", nmFile)
  if (length(minStart) > 1) {
    stop(
      "More than one set of termination messages in this file. This is not currently supported.\n"
    )
  } else
    if (length(minStart) == 0 || length(minEnd) == 0) {
      termMsg <- NULL
    } else {
      termMsg <- nmFile[(minStart + 1):(minEnd - 1)]
      termMsg <- substring(termMsg, 2)
    }
  extFileName <-
    paste0(sub("\\.\\w*$", "", fileName_read), ".ext")
  if (!file.exists(extFileName)) {
    stop(paste(
      "Could not find the raw results file (",
      extFileName,
      ") for  ",
      fileName_read,
      ".\n"
    ))
  } else {
    extData <- read.table(extFileName, skip = 1, header = TRUE)
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
