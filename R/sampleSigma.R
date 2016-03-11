#' Sample from the multivariate normal distribution using the SIGMA variance-covariance matrix to generate new sets of simulated EPSILONs from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315").
#' @param n Number of samples required.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' the estimated NONMEM SIGMA variance-covariance matrix. This provides \code{n} sets of EPSILON estimates
#' suitable for simulation of new datasets.
#'
#'
#' @examples
#'  sigDist <- sampleSigma("run315", 5000)

sampleSigma <- function(nmRun, n) {

  nmOutput <- readNM(nmRun)

  sigmas <- getSigma(nmOutput)

  mu   <- rep(0, times=ncol(sigmas))
  as.data.frame(mvrnorm(n=n, mu, Sigma=sigmas))
}
