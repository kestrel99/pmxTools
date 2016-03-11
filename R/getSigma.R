#' Extract SIGMA and related matrices from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{readNM}}.
#' @param output A flag specifying the matrix or matrices to be output. Valid flag values are "est" (the default),
#'  "se", "rse", "cor", "cse", or "all".
#'
#' @return A symmetrical matrix, or a list of symmetrical matrices if \code{all} is specified.
#'
#' \code{est} returns the estimated SIGMA variance-covariance matrix.
#' \code{se} returns the standard errors for the estimated SIGMA variance-covariance matrix.
#' \code{rse} returns the relative standard errors for the estimated SIGMA variance-covariance matrix (se/est*100).
#' \code{cor} returns the correlation matrix matrix.
#' \code{cse} returns the standard errors for the correlation matrix.
#' \code{all} returns all available SIGMA matrices.
#'
#' @examples
#'  nmOutput <- readNM("run315.xml")
#'  sigmas   <- getSigma(nmOutput)
#'  sigmaRSEs <- getSigma(nmOutput, "rse")

getSigma <- function(x, output = "est") {
  processMatrix <- function(worklist) {
    nrows <- length(worklist)
    m <-
      matrix(
        nrow = nrows,
        ncol = nrows,
        dimnames = list(
          paste("SIGMA", 1:nrows, sep = ""),
          paste("SIGMA", 1:nrows, sep = "")
        )
      )

    for (i in 1:nrows) {
      a <- as.numeric(unlist(worklist[i]))
      a <- a[1:length(a) - 1]
      a <- a[seq(1, length(a), by = 2)]
      m[i, 1:length(a)] <- a
    }

    for (i in 1:nrows) {
      for (j in 1:nrows) {
        if (is.na(m[i, j])) {
          m[i, j] <- m[j, i]
        }
      }
    }

    m[m == 1e+10] <- NA
    m
  }


  if (!(output %in% c("est", "se", "rse", "cor", "cse", "all"))) {
    stop("Please select a valid output option (est, se, rse, cor, cse, all).")
  }

  if (length(grep("row", names(x$nonmem$problem$estimation$sigma))) == 0) {
    ### single dimension
    if (output == "est") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$sigma)[1], dimnames=list("SIGMA1","SIGMA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "se") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$sigmase)[1], dimnames=list("SIGMA1","SIGMA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "cor") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$sigmac)[1], dimnames=list("SIGMA1","SIGMA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "cse") {
      o <- matrix(as.numeric(x$nonmem$problem$estimation$sigmacse)[1], dimnames=list("SIGMA1","SIGMA1"))
      o[o == 1e+10] <- NA
    }

    if (output == "rse") {
      m1 <- as.numeric(x$nonmem$problem$estimation$sigma)[1]
      m2 <- as.numeric(x$nonmem$problem$estimation$sigmase)[1]

      if(m2==1e+10) m2 <- NA

      m  <- matrix(m2 / m1 * 100, dimnames=list("SIGMA1","SIGMA1"))
      m[is.infinite(m)] <- NA
      o <- m
      o[o == 1e+10] <- NA
    }

    if (output == "all") {
      out <- list()

      out$Sigma <- matrix(as.numeric(x$nonmem$problem$estimation$sigma)[1], dimnames=list("SIGMA1","SIGMA1"))

      m <- as.numeric(x$nonmem$problem$estimation$sigmase)[1]
      if(m==1e+10) m <- NA
      out$SigmaSE <- matrix(m, dimnames=list("SIGMA1","SIGMA1"))

      m1 <- as.numeric(x$nonmem$problem$estimation$sigma)[1]
      m2 <- as.numeric(x$nonmem$problem$estimation$sigmase)[1]
      if(m2==1e+10) m2 <- NA
      m  <- matrix(m2 / m1 * 100, dimnames=list("SIGMA1","SIGMA1"))
      m[is.infinite(m)] <- NA
      out$SigmaRSE           <- m

      m <- as.numeric(x$nonmem$problem$estimation$sigmac)[1]
      if(m==1e+10) m <- NA
      out$SigmaCorrelation   <- matrix(m, dimnames=list("SIGMA1","SIGMA1"))

      m <- as.numeric(x$nonmem$problem$estimation$sigmacse)[1]
      if(m==1e+10) m <- NA
      out$SigmaCorrelationSE <- matrix(m, dimnames=list("SIGMA1","SIGMA1"))

      o <- out
    }

  } else {  ### matrix
    if (output == "est") {
      o <- processMatrix(x$nonmem$problem$estimation$sigma)
    }

    if (output == "se") {
      o <- processMatrix(x$nonmem$problem$estimation$sigmase)
    }

    if (output == "cor") {
      o <- processMatrix(x$nonmem$problem$estimation$sigmac)
    }

    if (output == "cse") {
      o <- processMatrix(x$nonmem$problem$estimation$sigmacse)
    }

    if (output == "rse") {
      m1 <- processMatrix(x$nonmem$problem$estimation$sigma)
      m2 <- processMatrix(x$nonmem$problem$estimation$sigmase)
      m  <- m2 / m * 100
      m[is.infinite(m)] <- NA
      o <- m
    }

    if (output == "all") {
      out <- list()
      out$Sigma              <-
        processMatrix(x$nonmem$problem$estimation$sigma)
      out$SigmaSE            <-
        processMatrix(x$nonmem$problem$estimation$sigmase)

      m1 <- processMatrix(x$nonmem$problem$estimation$sigma)
      m2 <- processMatrix(x$nonmem$problem$estimation$sigmase)
      m  <- m2 / m * 100
      m[is.infinite(m)] <- NA
      out$SigmaRSE           <- m

      out$SigmaCorrelation   <-
        processMatrix(x$nonmem$problem$estimation$sigmac)
      out$SigmaCorrelationSE <-
        processMatrix(x$nonmem$problem$estimation$sigmacse)

      o <- out
    }
  }
  o
}
