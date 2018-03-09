#' Sample from the multivariate normal distribution using the SIGMA variance-covariance matrix to generate new sets of simulated EPSILONs from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315").
#' @param n Number of samples required.
#' @param seed Random seed.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' the estimated NONMEM SIGMA variance-covariance matrix. This provides \code{n} sets of EPSILON estimates
#' suitable for simulation of new datasets.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  sigDist <- sample_sigma("run315", 5000, seed=740727)
#' }
#'
#' @export

sample_sigma <- function(nmRun, n, seed) {

  set.seed(seed)

  nmOutput <- read_nm(nmRun)

  sigmas <- get_sigma(nmOutput)

  mu   <- rep(0, times=ncol(sigmas))
  as.data.frame(MASS::mvrnorm(n=n, mu, Sigma=sigmas))
}
