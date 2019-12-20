#' Extract variability parameter estimates from a NONMEM output object.
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
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' \code{est} returns the estimated OMEGA variance-covariance matrix.
#' \code{se} returns the standard errors for the estimated OMEGA variance-covariance matrix.
#' \code{rse} returns the relative standard errors for the estimated OMEGA variance-covariance matrix (se/est*100).
#' \code{cor} returns the correlation matrix matrix.
#' \code{cse} returns the standard errors for the correlation matrix.
#' \code{95ci} returns the asymptotic 95\% confidence intervals for the elements of the OMEGA variance-covariance
#' matrix (est +/- 1.96*se).
#' \code{all} returns all available OMEGA matrices.
#'
#' @examples
#' \dontrun{
#'  nmOutput  <- read_nm("run315.xml")
#'  omegas    <- get_omega(nmOutput)
#'  omegaRSEs <- get_omega(nmOutput, "rse")
#' }
#'
#' @export

get_omega <- function (x, output = "est", sigdig = 6, sep = "-", est.step=NULL) 
{
  processMatrix <- function(worklist) {
    nrows <- length(worklist)
    m <- matrix(nrow = nrows, ncol = nrows, dimnames = list(paste("OMEGA", 
                                                                  1:nrows, sep = ""), paste("OMEGA", 1:nrows, sep = "")))
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
  
  if (length(grep("row", names(x$nonmem$problem[[ind_est]]$omega))) == 
      0) {
    if (output == "est") {
      o <- matrix(signif(as.numeric(x$nonmem$problem[[ind_est]]$omega)[1], 
                         sigdig), dimnames = list("OMEGA1", "OMEGA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "se") {
      o <- matrix(signif(as.numeric(x$nonmem$problem[[ind_est]]$omegase)[1], 
                         sigdig), dimnames = list("OMEGA1", "OMEGA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "cor") {
      o <- matrix(signif(as.numeric(x$nonmem$problem[[ind_est]]$omegac)[1], 
                         sigdig), dimnames = list("OMEGA1", "OMEGA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "cse") {
      o <- matrix(signif(as.numeric(x$nonmem$problem[[ind_est]]$omegacse)[1], 
                         sigdig), dimnames = list("OMEGA1", "OMEGA1"))
      o[o == 1e+10] <- NA
    }
    if (output == "rse") {
      m1 <- as.numeric(x$nonmem$problem[[ind_est]]$omega)[1]
      m2 <- as.numeric(x$nonmem$problem[[ind_est]]$omegase)[1]
      if (m2 == 1e+10) 
        m2 <- NA
      m <- matrix(m2/m1 * 100, dimnames = list("OMEGA1", 
                                               "OMEGA1"))
      m[is.infinite(m)] <- NA
      o <- signif(m, sigdig)
      o[o == 1e+10] <- NA
    }
    if (output == "95ci") {
      m1 <- as.numeric(x$nonmem$problem[[ind_est]]$omega)[1]
      m2 <- as.numeric(x$nonmem$problem[[ind_est]]$omegase)[1]
      if (m2 == 1e+10) 
        m2 <- NA
      oup <- signif(m1 + 1.96 * m2, sigdig)
      olo <- signif(m1 - 1.96 * m2, sigdig)
      o <- matrix(paste(olo, oup, sep = sep), dimnames = list("OMEGA1", 
                                                              "OMEGA1"))
    }
    if (output == "all") {
      out <- list()
      out$Omega <- matrix(signif(as.numeric(x$nonmem$problem[[ind_est]]$omega)[1], 
                                 sigdig), dimnames = list("OMEGA1", "OMEGA1"))
      m <- as.numeric(x$nonmem$problem[[ind_est]]$omegase)[1]
      if (m == 1e+10) 
        m <- NA
      out$OmegaSE <- matrix(signif(m, sigdig), dimnames = list("OMEGA1", 
                                                               "OMEGA1"))
      m1 <- as.numeric(x$nonmem$problem[[ind_est]]$omega)[1]
      m2 <- as.numeric(x$nonmem$problem[[ind_est]]$omegase)[1]
      if (m2 == 1e+10) 
        m2 <- NA
      m <- matrix(m2/m1 * 100, dimnames = list("OMEGA1", 
                                               "OMEGA1"))
      m[is.infinite(m)] <- NA
      out$OmegaRSE <- signif(m, sigdig)
      m1 <- as.numeric(x$nonmem$problem[[ind_est]]$omega)[1]
      m2 <- as.numeric(x$nonmem$problem[[ind_est]]$omegase)[1]
      if (m2 == 1e+10) 
        m2 <- NA
      oup <- signif(m1 + 1.96 * m2, sigdig)
      olo <- signif(m1 - 1.96 * m2, sigdig)
      out$Omega95CI <- paste(olo, oup, sep = sep)
      m <- as.numeric(x$nonmem$problem[[ind_est]]$omegac)[1]
      if (m == 1e+10) 
        m <- NA
      out$OmegaCorrelation <- matrix(signif(m, sigdig), 
                                     dimnames = list("OMEGA1", "OMEGA1"))
      m <- as.numeric(x$nonmem$problem[[ind_est]]$omegacse)[1]
      if (m == 1e+10) 
        m <- NA
      out$OmegaCorrelationSE <- matrix(signif(m, sigdig), 
                                       dimnames = list("OMEGA1", "OMEGA1"))
      o <- out
    }
  }
  else {
    if (output == "est") {
      o <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omega), 
                  sigdig)
    }
    if (output == "se") {
      o <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omegase), 
                  sigdig)
    }
    if (output == "cor") {
      o <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omegac), 
                  sigdig)
    }
    if (output == "cse") {
      o <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omegacse), 
                  sigdig)
    }
    if (output == "rse") {
      m1 <- processMatrix(x$nonmem$problem[[ind_est]]$omega)
      m2 <- processMatrix(x$nonmem$problem[[ind_est]]$omegase)
      m <- m2/m1 * 100
      m[is.infinite(m)] <- NA
      o <- signif(m, sigdig)
    }
    if (output == "95ci") {
      m1 <- processMatrix(x$nonmem$problem[[ind_est]]$omega)
      m2 <- processMatrix(x$nonmem$problem[[ind_est]]$omegase)
      mup <- signif(m1 + 1.96 * m2, sigdig)
      mlo <- signif(m1 - 1.96 * m2, sigdig)
      m <- matrix(paste(mlo, mup, sep = sep), nrow = nrow(m1), 
                  ncol = ncol(m2), dimnames = dimnames(m1))
      o <- m
    }
    if (output == "all") {
      out <- list()
      out$Omega <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omega), 
                          sigdig)
      out$OmegaSE <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omegase), 
                            sigdig)
      m1 <- processMatrix(x$nonmem$problem[[ind_est]]$omega)
      m2 <- processMatrix(x$nonmem$problem[[ind_est]]$omegase)
      m <- m2/m1 * 100
      m[is.infinite(m)] <- NA
      out$OmegaRSE <- signif(m, sigdig)
      m1 <- processMatrix(x$nonmem$problem[[ind_est]]$omega)
      m2 <- processMatrix(x$nonmem$problem[[ind_est]]$omegase)
      mup <- signif(m1 + 1.96 * m2, sigdig)
      mlo <- signif(m1 - 1.96 * m2, sigdig)
      m <- matrix(paste(mlo, mup, sep = sep), nrow = nrow(m1), 
                  ncol = ncol(m2), dimnames = dimnames(m1))
      out$Omega95CI <- m
      out$OmegaCorrelation <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omegac), 
                                     sigdig)
      out$OmegaCorrelationSE <- signif(processMatrix(x$nonmem$problem[[ind_est]]$omegacse), 
                                       sigdig)
      o <- out
    }
  }
  o
}