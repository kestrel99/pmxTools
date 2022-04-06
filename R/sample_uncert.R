#' Sample from the multivariate normal distribution to generate new sets of parameters from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315.xml").
#' @param n Number of samples required.
#' @param seed Random seed.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' NONMEM typical parameter estimates the NONMEM variance-covariance matrix (from the *.cov file). This
#' provides \code{n} sets of parameter estimates sampled from the uncertainty distribution, suitable
#' for simulation under model uncertainty.
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  nmMatrix <- sample_uncert("run315.xml", 5000, seed=740727)
#' }
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stringr str_replace
sample_uncert <- function(nmRun, n, seed) {
  set.seed(seed)
  nmOutput <- read_nm(nmRun)
  thetas <- get_theta(nmOutput)
  omegas <- get_omega(nmOutput)
  sigmas <- get_sigma(nmOutput)
  omList <- c()
  for (i in 1:nrow(omegas)) {
    for (j in 1:i) {
      omList <- c(omList, omegas[i, j])
    }
  }
  siList <- c()
  for (i in 1:nrow(sigmas)) {
    for (j in 1:i) {
      siList <- c(siList, sigmas[i, j])
    }
  }
  mu <- as.numeric(c(thetas, siList, omList))
  vcov <- read_nmcov(gsub(pattern = ".xml", replacement = ".cov", x = nmRun))
  as.data.frame(MASS::mvrnorm(n = n, mu, Sigma = vcov))
}
