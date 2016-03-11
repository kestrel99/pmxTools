#' Sample from the multivariate normal distribution using the OMEGA variance-covariance matrix to generate new sets of simulated ETAs from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315").
#' @param n Number of samples required.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' the estimated NONMEM OMEGA variance-covariance matrix. This provides \code{n} sets of ETA estimates
#' suitable for simulation of new patients.
#'
#'
#' @examples
#'  omDist <- sampleOmega("run315", 5000)

sampleOmega <- function(nmRun, n) {

  nmOutput <- readNM(nmRun)

  omegas <- getOmega(nmOutput)

  mu   <- rep(0, times=ncol(omegas))
  as.data.frame(mvrnorm(n=n, mu, Sigma=omegas))
}
