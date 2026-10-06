#' Sample from the multivariate normal distribution using the OMEGA variance-covariance matrix to generate new sets of simulated ETAs from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315").
#' @param n Number of samples required.
#' @param seed Random seed.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' the estimated NONMEM OMEGA variance-covariance matrix. This provides \code{n} sets of ETA estimates
#' suitable for simulation of new patients.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/solutions/technologies/nonmem})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  omDist <- sample_omega("run315", 5000, seed=740727)
#' }
#'
#' @export
#' @importFrom MASS mvrnorm
sample_omega <- function(nmRun, n, seed) {

  set.seed(seed)

  nmOutput <- read_nm(nmRun)

  omegas <- get_omega(nmOutput)

  mu   <- rep(0, times=ncol(omegas))
  as.data.frame(MASS::mvrnorm(n=n, mu, Sigma=omegas))
}
