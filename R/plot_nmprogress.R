#' Plot NONMEM parameter estimation by iteration.
#'
#' \code{plot_nmprogress} returns a plot or set of plots showing the evolution of 
#' parameter estimates by iteration.
#'
#' @param fileName A NONMEM output file prefix, without extension (e.g. 'run315').
#' @param fileExt  The file extension for NONMEM output, set to '.lst' by default.
#' @param metric  What to show in the plot. Allowed options are 'est' (the actual 
#' estimate) or 'perc' (the percentage change in the estimated or OFV since estimation began). 
#' Default is 'perc'.
#' @param lineCol  Line color. Default is '#902C10'.
#' @param idlineCol  Identity line color (only used if 'perc' metric is selected). Default is black.
#'
#' @return A set of plots.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' plot_nmprogress("run315")
#' plot_nmprogress("run315", ".nmlst")
#' }
#'
#' @import utils ggplot2
#' @export

plot_nmprogress <- function(fileName, fileExt = ".lst", metric="perc",
                       lineCol = "#902C10", idlineCol="black") {
  fileName <- paste(fileName, fileExt, sep = "")
  nmFile <-
    scan(fileName,
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
    paste(sub("\\.\\w*$", "", fileName), ".ext", sep = "")
  if (!file.exists(extFileName)) {
    stop(paste(
      "Could not find the raw results file (",
      extFileName,
      ") for  ",
      fileName,
      ".\n"
    ))
  } else {
    extData <- read.table(extFileName, skip = 1, header = T)
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

  pdata <- reshape(extData, idvar = "ITERATION", direction="long", varying = list(names(extData)[2:(length(extData))]))
  pdata$Parameter <- rep(names(extData)[2:(length(extData))], each=nrow(extData))
  pdata <- pdata[pdata$ITERATION >=0,]
  pdata$Parameter <- ordered(pdata$Parameter, unique(pdata$Parameter))
  
  a1 <- split(pdata, pdata$Parameter)
  a2 <- lapply(a1, function(x) {
    x$percDelta <- (x$THETA1 - x$THETA1[1]) / x$THETA1[1] * 100
    x
  })
  a3 <- do.call(rbind, a2)
  pdata <- a3
  
  drawplot <- function(var1, var2, ylab) {
    ggplot(pdata, aes(!! sym(var1), !! sym(var2))) +
      geom_line(col=lineCol) +
      facet_wrap(~ Parameter, scales = "free_y") +
      scale_x_continuous("Iteration") +
      scale_y_continuous(ylab)
  }
  
  if(metric!="perc") {
    drawplot("ITERATION", "THETA1", "Estimated value")
  } else {
    drawplot("ITERATION", "percDelta", "Change from starting value (%)")
  }
}

