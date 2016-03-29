#' Create a table of model parameter estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{readNM}}.
#'
#' @return A named vector of NONMEM model parameter estimates.
#'
#' @examples
#'  nmOutput <- readNM("run315.xml")
#'  estTab   <- getEstTable(nmOutput)

getEstTable <- function(x, thetaLabels=c(), omegaLabels=c(),
                        sigmaLabels=c(), sigdig=3) {

  theta     <- getTheta(x, sigdig=3)
  thetarse  <- getTheta(x, output="rse", sigdig=3)
  theta95ci <- getTheta(x, output="95ci", sigdig=3)

  omegaMatrix     <- getOmega(x, sigdig=3)
  omegaRSEMatrix  <- getOmega(x, sigdig=3, output="rse")
  omega95CIMatrix <- getOmega(x, sigdig=3, output="95ci")

  omegaL     <- c()
  omegaRSEL  <- c()
  omega95CIL <- c()
  omegaLab   <- c()
  for(n in 1:nrow(omegaMatrix)) {
    omegaL     <- c(omegaL, omegaMatrix[n,1:n])
    omegaRSEL  <- c(omegaRSEL, omegaRSEMatrix[n,1:n])
    omega95CIL <- c(omega95CIL, omega95CIMatrix[n,1:n])
    omegaLab   <- c(omegaLab, paste("OM", n, ",", 1:n, sep=""))
  }

  sigmaMatrix     <- getSigma(x, sigdig=3)
  sigmaRSEMatrix  <- getSigma(x, sigdig=3, output="rse")
  sigma95CIMatrix <- getSigma(x, sigdig=3, output="95ci")

  sigmaL     <- c()
  sigmaRSEL  <- c()
  sigma95CIL <- c()
  sigmaLab   <- c()
  for(n in 1:nrow(sigmaMatrix)) {
    sigmaL     <- c(sigmaL, sigmaMatrix[n,1:n])
    sigmaRSEL  <- c(sigmaRSEL, sigmaRSEMatrix[n,1:n])
    sigma95CIL <- c(sigma95CIL, sigma95CIMatrix[n,1:n])
    sigmaLab   <- c(sigmaLab, paste("SI", n, ",", 1:n, sep=""))
  }

  labCol <- paste("THETA", 1:length(theta),sep="")

  if(length(thetaLabels)==length(theta)) {
    labCol <- thetaLabels
  }

  idxOmegaDiag <- c()
  for (n in 1:nrow(omegaMatrix)){
    idxOmegaDiag <- c(idxOmegaDiag, sum(1:n))
  }
  idxOmegaOffDiag <- setdiff(1:length(omegaL), idxOmegaDiag)

  if(length(omegaLabels)==nrow(omegaMatrix)) {
    omegaLab[idxOmegaDiag] <- omegaLabels
  } else {
    omegaLab[idxOmegaDiag] <- paste("OMEGA", 1:nrow(omegaMatrix), sep="")
  }

  labCol <- c(labCol, omegaLab)

  idxSigmaDiag <- c()
  for (n in 1:nrow(sigmaMatrix)){
    idxSigmaDiag <- c(idxSigmaDiag, sum(1:n))
  }
  idxSigmaOffDiag <- setdiff(1:length(sigmaL), idxSigmaDiag)

  if(length(sigmaLabels)==nrow(sigmaMatrix)) {
    sigmaLab[idxSigmaDiag] <- sigmaLabels
  } else {
    sigmaLab[idxSigmaDiag] <- paste("SIGMA", 1:nrow(sigmaMatrix), sep="")
  }

  labCol <- c(labCol, sigmaLab)

  estCol <- c(as.numeric(theta),
              as.numeric(omegaL),
              as.numeric(sigmaL))

  rseCol <- c(as.numeric(thetarse),
              as.numeric(omegaRSEL),
              as.numeric(sigmaRSEL))

  CI95Col <- c(as.character(theta95ci),
               as.character(omega95CIL),
               as.character(sigma95CIL))

  shrEta <- rep("-", times=length(omegaL))
  shrEta[idxOmegaDiag] <- as.numeric(getShrinkage(x))

  shrCol <- c(rep("-", times=length(theta)),
              shrEta,
              as.numeric(getShrinkage(x, output="epsilon")))

  out <- data.frame(Parameter = labCol,
                    Estimate  = estCol,
                    RSE       = rseCol,
                    CI95      = CI95Col,
                    Shrinkage = shrCol)



  out
}
