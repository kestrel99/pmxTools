#' Extract residual variability parameter estimates from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param output A flag specifying the matrix or matrices to be output. Valid flag values are \code{est} (the default),
#'  \code{se}, \code{rse}, \code{cor}, \code{cse}, \code{95ci}, or \code{all}.
#' @param sigdig Specifies the number of significant digits to be provided (default=6).
#' @param sep Specifies the separator character to use for 95\% confidence intervals (default="-").
#' @param est.step Specifies which estimation step to return parameters from (default is the last).
#' 
#' @return A symmetrical matrix, or a list of symmetrical matrices if \code{all} is specified.
#'
#' \code{est} returns the estimated SIGMA variance-covariance matrix.
#' \code{se} returns the standard errors for the estimated SIGMA variance-covariance matrix.
#' \code{rse} returns the relative standard errors for the estimated SIGMA variance-covariance matrix (se/est*100).
#' \code{cor} returns the correlation matrix matrix.
#' \code{cse} returns the standard errors for the correlation matrix.
#' \code{95ci} returns the asymptotic 95\% confidence intervals for the elements of the SIGMA variance-covariance
#' matrix (est +/- 1.96*se).
#' \code{all} returns all available SIGMA matrices.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  nmOutput <- read_nm("run315.xml")
#'  sigmas   <- get_sigma(nmOutput)
#'  sigmaRSEs <- get_sigma(nmOutput, "rse")
#' }
#'
#' @export

get_sigma <- function (x, output = "est", sigdig = 6, sep = "-", est.step=NULL) 
{
  processMatrix <- function(worklist) {
    nrows <- length(worklist)
    m <- matrix(nrow = nrows, ncol = nrows, dimnames = list(paste("SIGMA", 
                                                                  1:nrows, sep = ""), paste("SIGMA", 1:nrows, sep = "")))
    for (i in 1:nrows) {
      a <- as.numeric(unlist(worklist[i]))
      # a <- a[1:length(a) - 1]
      # a <- a[seq(1, length(a), by = 2)]
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
  if (!(output %in% c("est", "se", "rse", "cor", "cse", "95ci", 
                      "all"))) {
    stop("Please select a valid output option (est, se, rse, 95ci, cor, cse, all).")
  }
  
  if(is.null(est.step)) {
    no_steps <- sum(stringr::str_count(names(x$nonmem$problem), "estimation"))
  } else {
    no_steps <- est.step
  }
  
  ind_est  <- match("estimation", names(x$nonmem$problem))-1+no_steps
  
  if (length(grep("row", names(x$nonmem$problem[[ind_est]]$sigma))) == 
      0) {
    if (output == "est") {
      o <- matrix(signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigma)), 
                         sigdig), dimnames = list("SIGMA1", "SIGMA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "se") {
      o <- matrix(signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmase)), 
                         sigdig), dimnames = list("SIGMA1", "SIGMA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "cor") {
      o <- matrix(signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmac)), 
                         sigdig), dimnames = list("SIGMA1", "SIGMA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "cse") {
      o <- matrix(signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmacse)), 
                         sigdig), dimnames = list("SIGMA1", "SIGMA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "rse") {
      m1 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigma))
      m2 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmase))
      if (m2 == 1e+10) 
        m2 <- NA
      m <- matrix(m2/m1 * 100, dimnames = list("SIGMA1", 
                                               "SIGMA1"))
      m[is.infinite(m)] <- NA
      o <- signif(m, sigdig)
      o[o == 1e+10] <- NA
    }
    if (output == "95ci") {
      m1 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigma))
      m2 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmase))
      if (m2 == 1e+10) 
        m2 <- NA
      oup <- signif(m1 + 1.96 * m2, sigdig)
      olo <- signif(m1 - 1.96 * m2, sigdig)
      o <- matrix(paste(olo, oup, sep = sep), dimnames = list("SIGMA1", 
                                                              "SIGMA1"))
    }
    if (output == "all") {
      out <- list()
      out$Sigma <- matrix(signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigma)), 
                                 sigdig), dimnames = list("SIGMA1", "SIGMA1"))
      m <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmase))
      if (m == 1e+10) 
        m <- NA
      out$SigmaSE <- matrix(signif(m, sigdig), dimnames = list("SIGMA1", 
                                                               "SIGMA1"))
      m1 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigma))
      m2 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmase))
      if (m2 == 1e+10) 
        m2 <- NA
      m <- matrix(m2/m1 * 100, dimnames = list("SIGMA1", 
                                               "SIGMA1"))
      m[is.infinite(m)] <- NA
      out$SigmaRSE <- signif(m, sigdig)
      m1 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigma))
      m2 <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmase))
      if (m2 == 1e+10) 
        m2 <- NA
      oup <- signif(m1 + 1.96 * m2, sigdig)
      olo <- signif(m1 - 1.96 * m2, sigdig)
      out$Sigma95CI <- paste(olo, oup, sep = sep)
      m <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmac))
      if (m == 1e+10) 
        m <- NA
      out$SigmaCorrelation <- matrix(signif(m, sigdig), 
                                     dimnames = list("SIGMA1", "SIGMA1"))
      m <- as.numeric(unlist(x$nonmem$problem[[ind_est]]$sigmacse))
      if (m == 1e+10) 
        m <- NA
      out$SigmaCorrelationSE <- matrix(signif(m, sigdig), 
                                       dimnames = list("SIGMA1", "SIGMA1"))
      o <- out
    }
  }
  else {
    if (output == "est") {
      o <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigma, unlist)), 
                  sigdig)
    }
    if (output == "se") {
      o <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmase, unlist)), 
                  sigdig)
    }
    if (output == "cor") {
      o <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmac, unlist)), 
                  sigdig)
    }
    if (output == "cse") {
      o <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmacse, unlist)), 
                  sigdig)
    }
    if (output == "rse") {
      m1 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigma, unlist))
      m2 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmase, unlist))
      m <- m2/m1 * 100
      m[is.infinite(m)] <- NA
      o <- signif(m, sigdig)
    }
    if (output == "95ci") {
      m1 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigma, unlist))
      m2 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmase, unlist))
      mup <- fmt_signif(m1 + 1.96 * m2, sigdig)
      mlo <- fmt_signif(m1 - 1.96 * m2, sigdig)
      m <- matrix(paste(mlo, mup, sep = sep), nrow = nrow(m1), 
                  ncol = ncol(m2), dimnames = dimnames(m1))
      o <- m
    }
    if (output == "all") {
      out <- list()
      out$Sigma <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigma, unlist)), 
                          sigdig)
      out$SigmaSE <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmase, unlist)), 
                            sigdig)
      m1 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigma, unlist))
      m2 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmase, unlist))
      m <- m2/m1 * 100
      m[is.infinite(m)] <- NA
      out$SigmaRSE <- signif(m, sigdig)
      m1 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigma, unlist))
      m2 <- processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmase, unlist))
      mup <- fmt_signif(m1 + 1.96 * m2, sigdig)
      mlo <- fmt_signif(m1 - 1.96 * m2, sigdig)
      m <- matrix(paste(mlo, mup, sep = sep), nrow = nrow(m1), 
                  ncol = ncol(m2), dimnames = dimnames(m1))
      out$Sigma95CI <- m
      out$SigmaCorrelation <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmac, unlist)), 
                                     sigdig)
      out$SigmaCorrelationSE <- signif(processMatrix(sapply(x$nonmem$problem[[ind_est]]$sigmacse, unlist)), 
                                       sigdig)
      o <- out
    }
  }
  o
}