#' Extract structural model parameter estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{readNM}}.
#' @param output A flag specifying the matrix or matrices to be output. Valid flag values are "est" (the default),
#'  "se", "rse", "cor", "cse", or "all".
#'
#' @return A symmetrical matrix, or a list of symmetrical matrices if \code{all} is specified.
#'
#' \code{est} returns the estimated OMEGA variance-covariance matrix.
#' \code{se} returns the standard errors for the estimated OMEGA variance-covariance matrix.
#' \code{rse} returns the relative standard errors for the estimated OMEGA variance-covariance matrix (se/est*100).
#' \code{cor} returns the correlation matrix matrix.
#' \code{cse} returns the standard errors for the correlation matrix.
#' \code{all} returns all available OMEGA matrices.
#'
#' @examples
#'  nmOutput <- readNM("run315.xml")
#'  omegas   <- getOmega(nmOutput)
#'  omegaRSEs <- getOmega(nmOutput, "rse")

getOmega <- function(x, output = "est") {

  processMatrix <- function(worklist) {
    nrows <- length(worklist)
    m <- matrix(nrow=4, ncol=4, dimnames=list(paste("OMEGA", 1:nrows, sep=""), paste("OMEGA", 1:nrows, sep="")))

    for(i in 1:nrows) {
      a <- as.numeric(unlist(worklist[i]))
      a <- a[1:length(a)-1]
      a <- a[seq(1, length(a), by=2)]
      m[i,1:length(a)] <- a
    }

    for(i in 1:nrows) {
      for(j in 1:nrows) {
        if(is.na(m[i,j])) {
          m[i,j] <- m[j,i]
        }
      }
    }

    m[m==1e+10] <- NA
    m
  }


  if(!(output %in% c("est","se","rse","cor","cse","all"))) {
    stop("Please select a valid output option (est, se, rse, cor, cse, all).")
  }

  if (length(grep("row", names(x$nonmem$problem$estimation$omega))) == 0) {
    ### single dimension
    if (output == "est") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$omega)[1], dimnames=list("OMEGA1","OMEGA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "se") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$omegase)[1], dimnames=list("OMEGA1","OMEGA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "cor") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$omegac)[1], dimnames=list("OMEGA1","OMEGA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "cse") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$omegacse)[1], dimnames=list("OMEGA1","OMEGA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "rse") {
      m1 <- as.numeric(x$nonmem$problem$estimation$omega)[1]
      m2 <- as.numeric(x$nonmem$problem$estimation$omegase)[1]

      if(m2==1e+10) m2 <- NA

      m  <- matrix(m2 / m1 * 100, dimnames=list("OMEGA1","OMEGA1"))
      m[is.infinite(m)] <- NA
      o <- m
      o[o == 1e+10] <- NA
    }

    if (output == "all") {
      out <- list()

      out$Omega <- matrix(as.numeric(x$nonmem$problem$estimation$omega)[1], dimnames=list("OMEGA1","OMEGA1"))

      m <- as.numeric(x$nonmem$problem$estimation$omegase)[1]
      if(m==1e+10) m <- NA
      out$OmegaSE <- matrix(m, dimnames=list("OMEGA1","OMEGA1"))

      m1 <- as.numeric(x$nonmem$problem$estimation$omega)[1]
      m2 <- as.numeric(x$nonmem$problem$estimation$omegase)[1]
      if(m2==1e+10) m2 <- NA
      m  <- matrix(m2 / m1 * 100, dimnames=list("OMEGA1","OMEGA1"))
      m[is.infinite(m)] <- NA
      out$OmegaRSE           <- m

      m <- as.numeric(x$nonmem$problem$estimation$omegac)[1]
      if(m==1e+10) m <- NA
      out$OmegaCorrelation   <- matrix(m, dimnames=list("OMEGA1","OMEGA1"))

      m <- as.numeric(x$nonmem$problem$estimation$omegacse)[1]
      if(m==1e+10) m <- NA
      out$OmegaCorrelationSE <- matrix(m, dimnames=list("OMEGA1","OMEGA1"))

      o <- out
    }

  } else {  ### matrix

    if(output=="est") {
      o <- processMatrix(x$nonmem$problem$estimation$omega)
    }

    if(output=="se") {
      o <- processMatrix(x$nonmem$problem$estimation$omegase)
    }

    if(output=="cor") {
      o <- processMatrix(x$nonmem$problem$estimation$omegac)
    }

    if(output=="cse") {
      o <- processMatrix(x$nonmem$problem$estimation$omegacse)
    }

    if(output=="rse") {
      m1 <- processMatrix(x$nonmem$problem$estimation$omega)
      m2 <- processMatrix(x$nonmem$problem$estimation$omegase)
      m  <- m2/m*100
      m[is.infinite(m)] <- NA
      o <- m
    }

    if(output=="all") {
      out <- list()
      out$Omega              <- processMatrix(x$nonmem$problem$estimation$omega)
      out$OmegaSE            <- processMatrix(x$nonmem$problem$estimation$omegase)

      m1 <- processMatrix(x$nonmem$problem$estimation$omega)
      m2 <- processMatrix(x$nonmem$problem$estimation$omegase)
      m  <- m2/m*100
      m[is.infinite(m)] <- NA
      out$OmegaRSE           <- m

      out$OmegaCorrelation   <- processMatrix(x$nonmem$problem$estimation$omegac)
      out$OmegaCorrelationSE <- processMatrix(x$nonmem$problem$estimation$omegacse)

      o <- out
    }}
  o
}
