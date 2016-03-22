#' Sample from the multivariate normal distribution to generate new sets of parameters from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315").
#' @param n Number of samples required.
#' @param seed Random seed.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' NONMEM typical parameter estimates the NONMEM variance-covariance matrix (from the *.cov file). This
#' provides \code{n} sets of parameter estimates sampled from the uncertainty distribution, suitable
#' for simulation under model uncertainty.
#'
#' @examples
#'  nmMatrix <- sampleUncert("run315", 5000, seed=740727)

sampleUncert <- function(nmRun, n, seed) {

  set.seed(seed)

  nmOutput <- readNM(nmRun)

  thetas <- getTheta(nmOutput)
  omegas <- getOmega(nmOutput)
  sigmas <- getSigma(nmOutput)

  omList <- c()
  for(i in 1:nrow(omegas)) {
    for(j in 1:i) {
      omList <- c(omList, omegas[i,j])
    }
  }

  siList <- c()
  for(i in 1:nrow(sigmas)) {
    for(j in 1:i) {
      siList <- c(siList, sigmas[i,j])
    }
  }

  mu   <- as.numeric(c(thetas, siList, omList))
  vcov <- readNMcov(nmRun)

  as.data.frame(mvrnorm(n=n, mu, Sigma=vcov))
}
