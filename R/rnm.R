#' Read NONMEM 7.2+ output into an R object.
#'
#' @param index The NONMEM model index, i.e. the numeric part of the filename assuming it follows the convention 'run123.mod'.
#' @param prefix The NONMEM model prefix, assuming it follows the convention 'run123.mod'. The default is \code{"run"}.
#' @param pathNM The path to the NONMEM output. This should not contain a trailing slash.
#' @param ndig Number of significant digits to use. The default is 3.
#' @param ndigB Number of significant digits to use. The default is 3.
#' @param ndigP Number of digits after the decimal point to use for percentages. The default is 1.
#' @param Pci Asymptotic confidence interval to apply when reporting parameter uncertainty. The default is 95.
#' @param ext NONMEM output file extension. The default is \code{".lst"}.
#' @param extmod NONMEM control stream file extension. The default is \code{"mod"}.
#' @param Pvalues Report P-values for parameters? The default is \code{TRUE}.
#' @param RawCI Report confidence intervals without estimate? The default is \code{FALSE}.
#' @param ... Additional arguments.
#'
#' @return A list containing information extracted from the NONMEM output.
#' @details The output list is composed of the following objects:
#' \itemize{
#'  \item{"Theta"}{A data frame describing the structural (fixed-effect) parameters, containing parameter name, estimated value, standard error (SE), coefficient of variation (CV), lower and upper confidence limits (CIL and CIU, based on \code{Pci}), and P-value, calculated as \code{2 * (1 - pnorm(abs(theta/theta.se)))}.}
#'  \item{"Eta"}{A data frame describing the interindividual random-effects parameters, containing estimated value, standard error (SE), coefficient of variation (CV, calculated as \code{abs(100*(SE/OMEGA))}), coefficient of variation (EtaCV, calculated as \code{100*sqrt(OMEGA)}), and shrinkage.}
#'  \item{"Epsilon"}{A data frame describing the residual random-effects parameters, containing estimated value, standard error (SE), coefficient of variation (CV, calculated as \code{abs(100*(SE/OMEGA))}), coefficient of variation (EtaCV, calculated as \code{100*sqrt(SIGMA)}), and shrinkage.}
#'  \item{"CorTheta"}{A data frame containing the correlation matrix for fixed effects (\code{"THETA"}).}
#'  \item{"CorOmega"}{A data frame containing the correlation matrix for interindividual random effects (\code{"OMEGA"}).}
#'  \item{"CorSigma"}{A data frame containing the correlation matrix for residual random effects (\code{"OMEGA"}).}
#'  \item{"OmegaMatrix"}{A data frame containing the \code{"OMEGA"} matrix.}
#'  \item{"SigmaMatrix"}{A data frame containing the \code{"OMEGA"} matrix.}
#'  \item{"CovMatrixTheta"}{A data frame containing the variance-covariance matrix for structural parameters (\code{THETA}).}
#'  \item{"CovMatrix"}{A data frame containing the complete variance-covariance matrix.}
#'  \item{"OFV"}{The objective function value.}
#'  \item{"ThetaString"}{A data frame containing all relevant fixed-effects parameter information, suitable for use in a table of parameter estimates. Contains parameter name, estimate, standard error, coefficient of variation, combined estimate and asymptotic confidence interval, and P-value.}
#'  \item{"EtaString"}{A data frame containing all relevant interindivudal random-effects parameter information, suitable for use in a table of parameter estimates. Contains parameter name, estimate (variance), standard error, coefficient of variation, percentage value (calculated as \code{100*sqrt(OMEGA)}), and shrinkage.}
#'  \item{"EpsString"}{A data frame containing all relevant residual random-effects parameter information, suitable for use in a table of parameter estimates. Contains parameter name, estimate (variance), standard error, coefficient of variation, percentage value (calculated as \code{100*sqrt(SIGMA)}), and shrinkage.}
#'  \item{"RunTime"}{Run time.}
#'  \item{"ConditionN"}{Condition number.}
#'  
#'  
#' }
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' 
#' @author Rik Schoemaker, \email{rik.schoemaker@@occams.com}
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' nmOutput <- rnm("run315.lst")
#' }
#'
#' @export


rnm <-
  function(index,
           prefix = "run",
           pathNM,
           ndig = 3,
           ndigB = 3,
           ndigP = 1,
           Pci = 95,
           ext = ".lst",
           extmod = ".mod",
           Pvalues = TRUE,
           RawCI = FALSE,
           ...) {
    
    
  FoldChange <- FALSE
  fchdig <- 2
  
  string2num <- function(x)
  {
    oldopts <- options(warn = -1)
    on.exit(options(oldopts))
    nc <- nchar(x)
    tmp <- substring(x, 1:nc, 1:nc)
    spc <- tmp == " "
    st <- !spc & c(T, spc[-nc])
    end <- !spc & c(spc[-1], T)
    as.numeric(substring(x, (1:nc)[st], (1:nc)[end]))
  }
  
  isReadable  <- function(filename)
  {
    if(!is.null(version$language)) {
      return(file.exists(filename)[1])
    }
  }
  
  FF <- function(x, dig) {
    formatC(x, digits = dig, format = "f")
  }
  
  FFP <-
    function(number, sigdig = 3) {
      res <-
        formatC(
          signif(number, digits = sigdig),
          digits = sigdig,
          format = "fg",
          flag = "#"
        )
      res <- sub('\\.$', "", res)
      res
    }
  
  ## from old xpose4 read.lst code
  readOutput <- function(filename) {
    filename <- file.path(pathNM, paste(prefix, index, ext, sep = ""))
    listfile <-
      scan(
        filename,
        sep = "\n",
        what = character(),
        quiet = TRUE,
        quote = ""
      )
    
    ## Find termination messages
    minimStart <- grep("#TERM:", listfile)
    minimEnd   <- grep("#TERE:", listfile)
    
    if (length(minimStart) > 1) {
      warning(
        "There seems to be more than one estimation step in this run. Only the last has been extracted.\n",
        call. = F
      )
      minimStart <- max(grep("#TERM:", listfile))
      minimEnd   <- max(grep("#TERE:", listfile))
    }
    
    if (length(minimStart) == 0 || length(minimEnd) == 0) {
      termes <- NULL
    } else {
      termes <- listfile[(minimStart + 1):(minimEnd - 1)]
      termes <- substring(termes, 2)
    }
    
    ## Figure out the name of the raw results file
    
    rawfile <- paste(sub("\\.\\w*$", '', filename), ".ext", sep = "")
    if (!isReadable(rawfile)) {
      stop(
        paste0(
          "Could not find the raw results file (",
          rawfile,
          ") for output file: ",
          filename,
          "\n"
        )
      )
    } else {
      extfile <-
        scan(
          rawfile,
          sep = "\n",
          what = character(),
          quiet = TRUE,
          quote = ""
        )
      TableLines <- max(grep("TABLE NO.", extfile))
      
      rawres <- read.table(rawfile, skip = TableLines, header = T)
    }
    
    ## Extract OFV
    ofv <- rawres$OBJ[rawres$ITERATION == -1000000000]
    
    ## Extract parameter estimates
    
    ## Get lines with relevant info
    finalEstLine <- as.numeric(rawres[rawres$ITERATION == -1000000000, ])
    
    ## Extract theta estimates
    thetas       <- finalEstLine[grep("THETA", names(rawres))]
    
    ## Extract omega estimates
    omindx <- grep("OMEGA", names(rawres))
    omega  <- list()
    seenOM <- NULL
    seenOM[1:100] <- 0
    
    for (om in omindx) {
      omnam <- names(rawres)[om]
      omcol <- as.numeric(sub("OMEGA\\.\\w*\\.", '', omnam, perl = TRUE))
      tmp1  <- sub("OMEGA\\.", '', omnam, perl = TRUE)
      omrow <- as.numeric(sub("\\.\\w*\\.$", '', tmp1, perl = TRUE))
      
      seenOM[omrow] <- seenOM[omrow] + 1
      if (seenOM[omrow] == 1) {
        omega[omrow][omcol] <- finalEstLine[om]
      } else {
        omega[[omrow]] <- c(omega[[omrow]], as.numeric(finalEstLine[om]))
        
      }
    }
    
    ## Extract sigma estimates
    siindx <- grep("SIGMA", names(rawres))
    sigma  <- list()
    seenSI <- NULL
    seenSI[1:100] <- 0
    
    if (length(siindx) != 0) {
      for (si in siindx) {
        sinam <- names(rawres)[si]
        sicol <- as.numeric(sub("SIGMA\\.\\w*\\.", '', sinam, perl = TRUE))
        tmp1  <- sub("SIGMA\\.", '', sinam, perl = TRUE)
        sirow <- as.numeric(sub("\\.\\w*\\.$", '', tmp1, perl = TRUE))
        
        seenSI[sirow] <- seenSI[sirow] + 1
        if (seenSI[sirow] == 1) {
          sigma[sirow][sicol] <- finalEstLine[si]
        } else {
          sigma[[sirow]] <- c(sigma[[sirow]], as.numeric(finalEstLine[si]))
          
        }
      }
    } else {
      sigma <- NULL
    }
    
    ## Extract condition number:
    c_line <- as.numeric(rawres[rawres$ITERATION == -1000000003, ])
    if (is.na(c_line[1])) {
      conditionN <- NULL
    } else {
      conditionN <- c_line[2]
    }
    
    ## Extract line with relevant info
    seEstLine <- as.numeric(rawres[rawres$ITERATION == -1000000001, ])
    
    if (is.na(seEstLine[1])) {
      sethetas <- NULL
      seomega <- NULL
      sesigma <- NULL
    } else {
      ## Extract theta estimates
      sethetas       <- seEstLine[grep("THETA", names(rawres))]
      for (it in 1:length(sethetas))
        if (sethetas[it] == 1.00000E+10)
          sethetas[it] <- 0
        
        ## Extract omega estimates
        omindx <- grep("OMEGA", names(rawres))
        seomega  <- list()
        seeseOM <- NULL
        seeseOM[1:100] <- 0
        
        for (om in omindx) {
          omnam <- names(rawres)[om]
          omcol <- as.numeric(sub("OMEGA\\.\\w*\\.", '', omnam, perl = TRUE))
          tmp1  <- sub("OMEGA\\.", '', omnam, perl = TRUE)
          omrow <- as.numeric(sub("\\.\\w*\\.$", '', tmp1, perl = TRUE))
          
          #RS modified: in specific cases we can arrive here and still have NaN omega estimates
          if (seEstLine[om] == 1.00000E+10 |
              is.na(seEstLine[om]))
            seEstLine[om] <- 0
          
          seeseOM[omrow] <- seeseOM[omrow] + 1
          if (seeseOM[omrow] == 1) {
            seomega[omrow][omcol] <- seEstLine[om]
          } else {
            seomega[[omrow]] <- c(seomega[[omrow]], as.numeric(seEstLine[om]))
            
          }
        }
        
        ## Extract sigma estimates
        siindx <- grep("SIGMA", names(rawres))
        sesigma  <- list()
        seenseSI <- NULL
        seenseSI[1:100] <- 0
        
        if (length(siindx) != 0) {
          for (si in siindx) {
            sinam <- names(rawres)[si]
            sicol <-
              as.numeric(sub("SIGMA\\.\\w*\\.", '', sinam, perl = TRUE))
            tmp1  <- sub("SIGMA\\.", '', sinam, perl = TRUE)
            sirow <- as.numeric(sub("\\.\\w*\\.$", '', tmp1, perl = TRUE))
            
            if (seEstLine[si] == 1.00000E+10)
              seEstLine[si] <- 0
            
            seenseSI[sirow] <- seenseSI[sirow] + 1
            if (seenseSI[sirow] == 1) {
              sesigma[sirow][sicol] <- seEstLine[si]
            } else {
              sesigma[[sirow]] <- c(sesigma[[sirow]], as.numeric(seEstLine[si]))
              
            }
          }
        } else {
          sesigma <- NULL
        }
        
    }
    
    ## Correlation matrix
    rawfile <- paste(sub("\\.\\w*$", '', filename), ".cor", sep = "")
    if (!isReadable(rawfile)) {
      warning(
        paste0(
          "Could not find the correlation matrix file (",
          rawfile,
          ") for list file: ",
          filename,
          "\n"
        )
      )
      cormat   <- NULL
      npcormat <- NULL
    } else {
      extfile <-
        scan(
          rawfile,
          sep = "\n",
          what = character(),
          quiet = TRUE,
          quote = ""
        )
      TableLines <- max(grep("TABLE NO.", extfile))
      
      rawres <- read.table(rawfile, skip = TableLines, header = T)
      cormat <-
        as.matrix.data.frame(rawres[1:length(thetas), 2:(length(thetas) + 1)])
      
      for (i in 1:length(thetas)) {
        cormat[i, i] <- 1
      }
      npcormat <- NULL
    }
    
    ## Covariance matrix
    rawfile <- paste(sub("\\.\\w*$", '', filename), ".cov", sep = "")
    if (!isReadable(rawfile)) {
      warning(
        paste0(
          "Could not find the covariance matrix file (",
          rawfile,
          ") for list file: ",
          filename,
          "\n"
        )
      )
      covmat  <- NULL
      fcovmat <- NULL
    } else {
      extfile <-
        scan(
          rawfile,
          sep = "\n",
          what = character(),
          quiet = TRUE,
          quote = ""
        )
      TableLines <- max(grep("TABLE NO.", extfile))
      rawres <- read.table(rawfile, skip = TableLines, header = T)
      covmat <- as.matrix.data.frame(rawres[1:length(thetas), 2:(length(thetas) + 1)])
      fcovmat <- as.matrix.data.frame(rawres[1:nrow(rawres), 2:(nrow(rawres)+1)])
    }
    
    ## Shrinkage
    rawfile <- paste(sub("\\.\\w*$", '', filename), ".shk", sep = "")
    if (!isReadable(rawfile)) {
      warning(paste0(
        "Could not find the shrinkage file (",
        rawfile,
        ") for list file: ",
        filename,
        "\n"
      ))
      shkmat <- NULL
    } else {
      extfile <-
        scan(
          rawfile,
          sep = "\n",
          what = character(),
          quiet = TRUE,
          quote = ""
        )
      TableLines <- max(grep("TABLE NO.", extfile))
      rawres <- read.table(rawfile, skip = TableLines, header = TRUE)
      shkmat <- as.matrix.data.frame(rawres)
      shkmat <- rawres
      
      shrink <- shkmat[shkmat$TYPE == 4, -c(1)]
      epsshrink <- shkmat[shkmat$TYPE == 5, -c(1, 2)]
      epsshrink <- epsshrink[epsshrink != 0]
      if (length(epsshrink) == 0) {
        epsshrink = NA
      }
      names(epsshrink) <- paste("EPS", seq(1, length(epsshrink)), sep = ".")
    }
    
    ret.list <- list(term = termes, ofv = ofv, 
                     thetas = thetas, omega = omega, sigma = sigma, 
                     sethetas = sethetas, seomegas = seomega, sesigmas = sesigma,
                     cormat = cormat, npcormat = npcormat, covmat = covmat, fcovmat = fcovmat,
                     shrink = shrink, epsshrink = epsshrink, conditionN = conditionN)
  }
  
  ncrit <- qnorm((1 - (1 - Pci / 100) / 2))
  ncrit90 <- qnorm((1 - (1 - 90 / 100) / 2))
  
  fres <-
    as.data.frame(scan(
      file = file.path(pathNM, paste(prefix, index, ext, sep = "")),
      sep = "\t",
      what = character(),
      quote = "", quiet=T
    ),
    stringsAsFactors = F)[, 1]
  
  fstart <- fres[grep("Run Start Time:", fres)]
  fstop  <- fres[grep("Run Stop Time:", fres)]

  if (length(fstop) > 0) {
    Tstart <- chron::chron(times. = substring(fstart, 17, 30),
                    format = c("h:m:s"))
    Tstop <- chron::chron(times. = substring(fstop, 17, 30), format = c("h:m:s"))
    Time <- (as.numeric(Tstop - Tstart) * 60 * 60 * 24)
    #if the run starts before and ends after midnight, this breaks down :)
    if (Time < 0) {
      Time <- Time + 60 * 60 * 24
    }
  } else {
    Time <- 0
  }
  
  if (exists("THETA"))
    rm(THETA, inherits = TRUE)
  if (exists("ETAD"))
    rm(ETAD, inherits = TRUE)
  if (exists("EPSD"))
    rm(EPSD, inherits = TRUE)
  
  f <- readOutput(filename = file.path(pathNM, paste(prefix, index, ext, sep = "")))  # from xpose4
  
  name <- paste(prefix, index, sep = "")
  fobj <-
    paste(
      "MINIMUM VALUE OF OBJECTIVE FUNCTION: " ,
      formatC(as.numeric(f$ofv), digits = 3, format = "f"),
      sep = ""
    )
  
  minimisation <- c(f$term, "", fobj)
  # Get number of thetas
  nth <- length(f$thetas)
  # Get number of etas
  net <- length(f$omega)
  # Get number of epsilons
  nep <- length(f$sigma)
  
  #Read names of THETAs
  filename <- file.path(pathNM, paste(prefix, index, ext, sep = ""))
  rawfile <- paste(sub("\\.\\w*$", '', filename), extmod, sep = "")
  if (!isReadable(rawfile)) {
    warning(
      paste0(
        "Could not find the control stream (",
        rawfile,
        ") for list file: ",
        filename,
        "\n"
      )
    )
    cormat   <- NULL
    npcormat <- NULL
  } else {
  NM <-
    as.data.frame(scan(
      file = file.path(pathNM, paste(name, extmod, sep = "")),
      sep = "\n",
      what = character(),
      quiet = T
    ), quote = "")
  names(NM)[1] <- "V1"
  NM$Parameter <-
    substr(NM$V1, ifelse(regexpr(";", NM$V1) == -1, 1000, regexpr(";", NM$V1) +
                           1), 1000)
  # lines with $THETA
  NM$TS1 <- as.numeric(regexpr("\\$THETA", NM$V1))
  f$Parameter <- trimws(NM[NM$TS1 == 1, "Parameter"])
  if (f$Parameter[1] == "")
    f$Parameter <- paste("THETA", seq(1:nth), sep = "")
  if (length(f$Parameter) != nth)
    f$Parameter <- paste("THETA", seq(1:nth), sep = "")
  }
  
  # omega as matrix:
  if (!is.null(f$omega)) {
    net <- length(f$omega)
    omegaM <- matrix(0, ncol = net, nrow = net)
    k <- 0
    for (i in 1:net) {
      k <- k + 1
      for (j in 1:k)
        omegaM[i, j] <- f$omega[[i]][j]
    }
    k <- 1
    if (net > 1) {
      for (i in 1:(net - 1)) {
        k <- k + 1
        for (j in k:net)
          omegaM[i, j] <- omegaM[j, i]
      }
    }
    omegaM <- as.matrix.data.frame(omegaM)
  } else {
    omegaM <- NULL
  }
  
  # sigma as matrix:
  if (!is.null(f$sigma)) {
    nep <- length(f$sigma)
    sigmaM <- matrix(0, ncol = nep, nrow = nep)
    k <- 0
    for (i in 1:nep) {
      k <- k + 1
      for (j in 1:k)
        sigmaM[i, j] <- f$sigma[[i]][j]
    }
    k <- 1
    if (nep > 1) {
      for (i in 1:(nep - 1)) {
        k <- k + 1
        for (j in k:nep)
          sigmaM[i, j] <- sigmaM[j, i]
      }
    }
    sigmaM <- as.matrix.data.frame(sigmaM)
  } else {
    sigmaM <- NULL
  }
  
  if (!is.null(f$sethetas)) {
    f$cvthetas     <- abs((f$sethetas / f$thetas) * 100)
    f$cithetasl     <- f$thetas - ncrit * f$sethetas
    f$cithetasu     <- f$thetas + ncrit * f$sethetas
    f$pvalue <- 2 * (1 - pnorm(abs(f$thetas / f$sethetas)))
    
    THETAs <-
      data.frame(
        stringsAsFactors = FALSE,
        Parameter = f$Parameter,
        Estimate = FFP(f$thetas, ndig),
        SE = FFP(f$sethetas, ndig),
        CV = paste(formatC(
          f$cvthetas, format = "f", digits = 1
        ), "%", sep = ""),
        CI = paste(FFP(f$cithetasl, ndig),
          "/",
          FFP(f$cithetasu, ndig),
          sep = ""
        ),
        CIE = paste(
          FFP(f$thetas, ndig),
          " (",
          FFP(f$cithetasl, ndig),
          "/",
          FFP(f$cithetasu, ndig),
          ")",
          sep = ""
        ),
        CIT = paste(
          FFP(10 ** f$thetas, ndigB),
          " (",
          FFP(10 ** f$cithetasl, ndigB),
          "/",
          FFP(10 ** f$cithetasu, ndigB),
          ")",
          sep = ""
        ),
        lCIT = paste(
          FFP(exp(f$thetas), ndigB),
          " (",
          FFP(exp(f$cithetasl), ndigB),
          "/",
          FFP(exp(f$cithetasu), ndigB),
          ")",
          sep = ""
        ),
        lCITp = paste(
          FF(100 * (exp(f$thetas) - 1), ndigP),
          "% (",
          FF(100 * (exp(f$cithetasl) - 1), ndigP),
          "%/",
          FF(100 * (exp(f$cithetasu) - 1), ndigP),
          "%)",
          sep = ""
        ),
        fchCss = paste(
          FF(1 / (exp(f$thetas)), fchdig),
          " (",
          FF(1 / (exp(
            f$thetas + ncrit90 * f$sethetas
          )), fchdig),
          "/",
          FF(1 / (exp(
            f$thetas - ncrit90 * f$sethetas
          )), fchdig),
          ")",
          sep = ""
        ),
        
        logitCIT = paste(FFP((
          exp(f$thetas) / (1 + exp(f$thetas))
        ), ndigB),
        " (", FFP((
          exp(f$cithetasl) / (1 + exp(f$cithetasl))
        ), ndigB),
        "/", FFP((
          exp(f$cithetasu) / (1 + exp(f$cithetasu))
        ), ndigB), ")", sep = ""),
        P = ifelse(
          is.na(f$pvalue),
          "",
          ifelse(f$pvalue < 0.00001, "< 0.00001", FF(f$pvalue, 5))
        )
      )
    
    THETAs$tCIT <- THETAs$CI
    THETAs$tCIT[regexpr("log", THETAs$Parameter) > -1]     <- THETAs$lCIT[regexpr("log", THETAs$Parameter) > -1]
    THETAs$tCIT[regexpr("logit", THETAs$Parameter) > -1]   <- THETAs$logitCIT[regexpr("logit", THETAs$Parameter) > -1]
    THETAs$tCIT[regexpr("10log", THETAs$Parameter) > -1]   <- THETAs$CIT[regexpr("10log", THETAs$Parameter) > -1]
    THETAs$tCIT[regexpr("logP", THETAs$Parameter) > -1]    <- THETAs$lCITp[regexpr("logP", THETAs$Parameter) > -1]
    THETAs$fchCss[regexpr("logP", THETAs$Parameter) <= -1] <- ""
    
    THETAs$CI <- THETAs$tCIT
    
    THETAs$tCITE <- THETAs$CIE
    THETAs$tCITE[regexpr("log", THETAs$Parameter) > -1]     <- THETAs$lCIT[regexpr("log", THETAs$Parameter) > -1]
    THETAs$tCITE[regexpr("logit", THETAs$Parameter) > -1]   <- THETAs$logitCIT[regexpr("logit", THETAs$Parameter) > -1]
    THETAs$tCITE[regexpr("10log", THETAs$Parameter) > -1]   <- THETAs$CIT[regexpr("10log", THETAs$Parameter) > -1]
    THETAs$tCITE[regexpr("logP", THETAs$Parameter) > -1]    <- THETAs$lCITp[regexpr("logP", THETAs$Parameter) > -1]
    THETAs$fchCss[regexpr("logP", THETAs$Parameter) <= -1] <- ""
    
    THETAs$CIE <- THETAs$tCITE
    
    THETAs$lCIT     <- NULL
    THETAs$lCITp    <- NULL
    THETAs$logitCIT <- NULL
    THETAs$tCIT     <- NULL
    THETAs$tCITE    <- NULL
    THETAs$CIT    <- NULL
    
    THETAs$SE[is.na(f$cvthetas)] <- "(Fixed)"
    THETAs$CV[is.na(f$cvthetas)] <- " "
    THETAs$CI[is.na(f$cvthetas)] <- " "
    THETAs$CIE[is.na(f$cvthetas)] <- " "
    names(THETAs)[5] <- paste(as.character(Pci), "% CI", sep = "")
    NameCI <- paste("Estimate", " (", names(THETAs)[5], ")", sep = "")
    names(THETAs)[6] <- NameCI
    if (RawCI == FALSE)   {
      THETAs[, 5] <- NULL
    }
    THETA  <-
      data.frame(
        Parameter = f$Parameter,
        Estimate = f$thetas,
        SE = f$sethetas,
        CV = f$cvthetas,
        CIL = f$cithetasl,
        CIU = f$cithetasu,
        P = f$pvalue
      )
    if (Pvalues == FALSE) {
      THETAs$P <- NULL
    }
    THETAs$"Fold change in Css (90% CI)" <- THETAs$fchCss
    THETAs$fchCss <- NULL
    if (FoldChange == FALSE) {
      THETAs$"Fold change in Css (90% CI)" <- NULL
    }
  }
  
  if (is.null(f$sethetas)) {
    empty <- as.vector(mode = "numeric", matrix(data = NA, nrow = nth, ncol = 1))
    emptyS <- as.vector(mode = "character", matrix(data = " ", nrow = nth, ncol = 1))
    
    THETA  <-
      data.frame(
        Parameter = f$Parameter,
        Estimate = f$thetas,
        SE = empty,
        CV = empty,
        CIL = empty,
        CIU = empty
      )
    
    THETAs <-
      data.frame(
        stringsAsFactors = FALSE,
        Parameter = f$Parameter,
        Estimate = FFP(f$thetas, ndig),
        SE = emptyS,
        CV = emptyS,
        CI = FFP(f$thetas, ndig),
        CIT = FFP(10 ** f$thetas, ndigB),
        lCIT = FFP(exp(f$thetas), ndigB),
        lCITp = paste(FF(100 * (exp(
          f$thetas
        ) - 1), ndigP), "%", sep = ""),
        fchCss = FF(1 / (exp(f$thetas)), fchdig),
        logitCIT = FFP((exp(f$thetas) / (1 + exp(
          f$thetas
        ))), ndigB),
        P = emptyS
      )
    
    THETAs$tCIT <- THETAs$CI
    THETAs$tCIT[regexpr("log", THETAs$Parameter) > -1]     <-
      THETAs$lCIT[regexpr("log", THETAs$Parameter) > -1]
    THETAs$tCIT[regexpr("logit", THETAs$Parameter) > -1]   <-
      THETAs$logitCIT[regexpr("logit", THETAs$Parameter) > -1]
    THETAs$tCIT[regexpr("10log", THETAs$Parameter) > -1]   <-
      THETAs$CIT[regexpr("10log", THETAs$Parameter) > -1]
    THETAs$tCIT[regexpr("logP", THETAs$Parameter) > -1]    <-
      THETAs$lCITp[regexpr("logP", THETAs$Parameter) > -1]
    THETAs$fchCss[regexpr("logP", THETAs$Parameter) <= -1] <- ""
    
    THETAs$CIT <- THETAs$tCIT
    
    THETAs$lCIT     <- NULL
    THETAs$lCITp    <- NULL
    THETAs$logitCIT <- NULL
    THETAs$tCIT     <- NULL
    THETAs$CI       <- " "
    
    names(THETAs)[5] <- paste(as.character(Pci), "% CI", sep = "")
    names(THETAs)[6] <- "Back-transformed estimate"
    if (RawCI == FALSE)   {
      THETAs[, 5] <- NULL
    }
    if (Pvalues == FALSE) {
      THETAs$P <- NULL
    }
    THETAs$"Fold change in Css (90% CI)" <- THETAs$fchCss
    THETAs$fchCss <- NULL
    if (FoldChange == FALSE) {
      THETAs$"Fold change in Css (90% CI)" <- NULL
    }
  }
  
  f$fetD   <-
    as.vector(mode = "numeric", matrix(data = NA, nrow = net, ncol = 1))
  f$fetSED <- f$fetD
  f$fetCVD <- f$fetD
  for (i in 1:net) {
    f$fetD[i] <- f$omega[[i]][i]
    f$fetDsqrt[i] <- 100 * sqrt(f$fetD[i])
  }
  if (!is.null(f$seomegas)) {
    for (i in 1:net) {
      f$fetSED[i] <- f$seomegas[[i]][i]
      if (f$fetD[i] == 0) {
        f$fetCVD[i] <- 0
      }
      else {
        f$fetCVD[i] <- abs(100 * (f$fetSED[i] / f$fetD[i]))
      }
    }
  }
  
  f$fepD   <-
    as.vector(mode = "numeric", matrix(data = NA, nrow = nep, ncol = 1))
  f$fepSED <- f$fepD
  f$fepCVD <- f$fepD
  for (i in 1:nep) {
    f$fepD[i] <- f$sigma[[i]][i]
    f$fepDsqrt[i] <- 100 * sqrt(f$fepD[i])
  }
  if (!is.null(f$sesigmas)) {
    for (i in 1:nep) {
      f$fepSED[i] <- f$sesigmas[[i]][i]
      if (f$fepD[i] == 0) {
        f$fepCVD[i] <- 0
      }
      else {
        f$fepCVD[i] <- abs(100 * (f$fepSED[i] / f$fepD[i]))
      }
    }
  }
  
  Shrinkage <- as.data.frame(t(f$shrink))
  SubPops <- dim(Shrinkage)[2]
  if (SubPops > 1) {
    nams <- paste("Shrink Pop", Shrinkage[1,], sep = "")
  } else {
    nams <- "Shrinkage"
  }
  Shrinkage <- Shrinkage[-1,]
  Shrinkage <- matrix(unlist(Shrinkage), ncol = SubPops)
  ShrinkageS <-
    paste(formatC(Shrinkage, format = "f", digits = 1), "%", sep = "")
  Shrinkage <- as.data.frame(matrix(Shrinkage, ncol = SubPops))
  names(Shrinkage) <- nams
  ShrinkageS <- as.data.frame(matrix(ShrinkageS, ncol = SubPops))
  names(ShrinkageS) <- nams
  
  EShrinkage <- as.data.frame(t(f$epsshrink))
  ESubPops <- dim(EShrinkage)[2]
  if (ESubPops > 1) {
    nams <- paste("Shrink Pop", EShrinkage[1,], sep = "")
  } else {
    nams <- "Shrinkage"
  }
  #EShrinkage <- EShrinkage[-1,]
  EShrinkage <- matrix(unlist(EShrinkage), ncol = ESubPops)
  EShrinkageS <-
    paste(formatC(EShrinkage, format = "f", digits = 1), "%", sep = "")
  EShrinkage <- as.data.frame(matrix(EShrinkage, ncol = ESubPops))
  names(EShrinkage) <- nams
  EShrinkageS <- as.data.frame(matrix(EShrinkageS, ncol = ESubPops))
  names(EShrinkageS) <- nams
  
  EpsShrinkage <- EShrinkage
  
  
  if (!is.null(f$seomegas)) {
    ETADs <-
      data.frame(
        stringsAsFactors = FALSE,
        Eta = 1:net,
        Estimate = FFP(f$fetD, ndig),
        SE = FFP(f$fetSED, ndig),
        CV = paste(formatC(
          f$fetCVD, format = "f", digits = 1
        ), "%", sep = "")
        ,
        EtaCV = paste(formatC(
          f$fetDsqrt, format = "f", digits = 1
        ), "%", sep = ""),
        Shrinkage = ShrinkageS
      )
    
    ETADs$SE[f$fetCVD <= 0] <- "(Fixed)"
    ETADs$CV[f$fetCVD <= 0] <- " "
    ETADs$EtaCV[f$fetCVD <= 0] <- " "
    
    ETAD  <-
      data.frame(
        Eta = 1:net,
        Estimate = f$fetD,
        SE = f$fetSED,
        CV = f$fetCVD,
        EtaCV = f$fetDsqrt,
        Shrinkage = Shrinkage
      )
  }
  
  if (is.null(f$seomegas)) {
    ETADs <-
      data.frame(
        Eta = 1:net,
        Estimate = FFP(f$fetD, ndig),
        SE = f$fetSED,
        CV = f$fetCVD,
        EtaCV = paste(formatC(
          f$fetDsqrt, format = "f", digits = 1
        ), "%", sep = ""),
        Shrinkage = ShrinkageS
      )
    ETAD  <-
      data.frame(
        Eta = 1:net,
        Estimate = f$fetD,
        SE = f$fetSED,
        CV = f$fetCVD,
        EtaCV = f$fetDsqrt,
        Shrinkage = Shrinkage
      )
  }
  
  names(ETADs)[5] <- '100*sqrt(Omega)'
  names(ETADs)[seq(6, dim(ETADs)[2])] <- nams
  
  if (!is.null(f$sesigmas)) {
    EPSDs <-
      data.frame(
        stringsAsFactors = FALSE,
        Epsilon = 1:nep,
        Estimate = FFP(f$fepD, ndig),
        SE = FFP(f$fepSED, ndig),
        CV = paste(formatC(
          f$fepCVD, format = "f", digits = 1
        ), "%", sep = "")
        ,
        EpsCV = paste(formatC(
          f$fepDsqrt, format = "f", digits = 1
        ), "%", sep = ""),
        Shrinkage = EShrinkageS
      )
    
    EPSDs$SE[f$fepCVD <= 0] <- "(Fixed)"
    EPSDs$CV[f$fepCVD <= 0] <- " "
    EPSDs$EpsCV[f$fepCVD <= 0] <- " "
    
    EPSD  <-
      data.frame(
        Epsilon = 1:nep,
        Estimate = f$fepD,
        SE = f$fepSED,
        CV = f$fepCVD,
        EpsCV = f$fepDsqrt,
        Shrinkage = EShrinkage
      )
  }
  
  if (is.null(f$sesigmas)) {
    EPSDs <-
      data.frame(
        Epsilon = 1:nep,
        Estimate = FFP(f$fepD, ndig),
        SE = f$fepSED,
        CV = f$fepCVD,
        EpsCV = paste(formatC(
          f$fepDsqrt, format = "f", digits = 1
        ), "%", sep = ""),
        Shrinkage = EShrinkageS
      )
    EPSD  <-
      data.frame(
        Epsilon = 1:nep,
        Estimate = f$fepD,
        SE = f$fepSED,
        CV = f$fepCVD,
        EpsCV = f$fepDsqrt,
        Shrinkage = EShrinkage
      )
  }
  
  names(EPSDs)[5] <- '100*sqrt(Omega)'
  names(EPSDs)[seq(6, dim(EPSDs)[2])] <- nams
  
  
  offdiag <- FALSE
  f$CorOmega   <- matrix(data = NA, nrow = net, ncol = net)
  for (i in 1:net) {
    for (j in 1:i) {
      VARi <- sqrt(f$omega[[i]][i])
      VARj <- sqrt(f$omega[[j]][j])
      if (VARi == 0 | VARj == 0) {
        f$CorOmega[i, j] <- 0
      }
      else {
        f$CorOmega[i, j] <- f$omega[[i]][j] / (VARi * VARj)
      }
      if ((i != j) & (f$CorOmega[i, j] != 0)) {
        offdiag <- TRUE
      }
    }
  }
  k <- 1
  if (net > 1) {
    for (i in 1:(net - 1)) {
      k <- k + 1
      for (j in k:net) {
        f$CorOmega[i, j] <- f$CorOmega[j, i]
      }
    }
  }
  tmp   <- as.data.frame(matrix(data = NA, nrow = net, ncol = net))
  for (i in 1:net) {
    for (j in 1:net) {
      tmp[i, j] <- formatC(f$CorOmega[i, j], format = "f", digits = 3)
    }
  }
  rownames(tmp) <- paste("Eta", 1:net, sep = "")
  colnames(tmp) <- paste("Eta", 1:net, sep = "")
  f$CorOmega <- tmp
  
  offdiag <- FALSE
  f$CorSigma   <- matrix(data = NA, nrow = nep, ncol = nep)
  for (i in 1:nep) {
    for (j in 1:i) {
      VARi <- sqrt(f$sigma[[i]][i])
      VARj <- sqrt(f$sigma[[j]][j])
      if (VARi == 0 | VARj == 0) {
        f$CorSigma[i, j] <- 0
      }
      else {
        f$CorSigma[i, j] <- f$sigma[[i]][j] / (VARi * VARj)
      }
      if ((i != j) & (f$CorSigma[i, j] != 0)) {
        offdiag <- TRUE
      }
    }
  }
  k <- 1
  if (nep > 1) {
    for (i in 1:(nep - 1)) {
      k <- k + 1
      for (j in k:nep) {
        f$CorSigma[i, j] <- f$CorSigma[j, i]
      }
    }
  }
  tmp   <- as.data.frame(matrix(data = NA, nrow = nep, ncol = nep))
  for (i in 1:nep) {
    for (j in 1:nep) {
      tmp[i, j] <- formatC(f$CorSigma[i, j], format = "f", digits = 3)
    }
  }
  rownames(tmp) <- paste("Eps", 1:nep, sep = "")
  colnames(tmp) <- paste("Eps", 1:nep, sep = "")
  f$CorSigma <- tmp
  
  params <-
    list(
      Theta = THETA,
      Eta = ETAD,
      Epsilon = EPSD,
      CorTheta = f$cormat,
      CorOmega = f$CorOmega,
      CorSigma = f$CorSigma,
      #NPCorMatrix = f$npcormat,
      OmegaMatrix = omegaM,
      SigmaMatrix = sigmaM,
      CovMatrixTheta = f$covmat,
      CovMatrix = f$fcovmat,
      OFV = f$ofv,
      ThetaString = THETAs,
      EtaString = ETADs,
      EpsString = EPSDs,
      #EpsShrink = f$epsshrink,
      RunTime = Time,
      conditionN = f$conditionN
    )
  return(params)
  
}
