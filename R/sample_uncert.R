#' Sample from the multivariate normal distribution to generate new sets of parameters from NONMEM output.
#'
#' @param nmRun Root filename for the NONMEM run (e.g. "run315.xml").
#' @param n Number of samples required.
#' @param seed Random seed.
#'
#' @return A data frame containing \code{n} samples from the multivariate normal distribution, using
#' NONMEM typical parameter estimates the NONMEM variance-covariance matrix (from the *.cov file). This
#' provides \code{n} sets of parameter estimates sampled from the uncertainty distribution, suitable
#' for simulation under model uncertainty. Columns are named after the parameters in the
#' \code{.cov} file (e.g. \code{THETA1}, \code{SIGMA(1,1)}, \code{OMEGA(2,1)}).
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/solutions/technologies/nonmem})
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
  # Name each estimate as NONMEM labels it in the .cov file, so that the means
  # can be matched to the covariance matrix by name rather than by position
  mu <- c(
    stats::setNames(as.numeric(thetas), paste0("THETA", seq_along(thetas))),
    lower_triangle(sigmas, "SIGMA"),
    lower_triangle(omegas, "OMEGA")
  )
  vcov <- read_nmcov(sub(pattern = "\\.xml$", replacement = ".cov", x = nmRun))
  labels <- rownames(vcov)
  if (!identical(labels, colnames(vcov)) || !setequal(labels, names(mu)) || length(labels) != length(mu)) {
    stop("The parameters in the covariance matrix do not match the estimates in ", nmRun, ".", call. = FALSE)
  }
  mu <- mu[labels]
  draws <- MASS::mvrnorm(n = n, mu = mu, Sigma = vcov)
  # mvrnorm() returns a vector for a single draw
  as.data.frame(matrix(draws, nrow = n, dimnames = list(NULL, labels)), optional = TRUE)
}

# Lower-triangle elements of an OMEGA or SIGMA matrix, row by row, named as in
# NONMEM's .cov file (e.g. "OMEGA(2,1)")
lower_triangle <- function(m, prefix) {
  idx <- which(lower.tri(m, diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
  stats::setNames(as.numeric(m[idx]), sprintf("%s(%d,%d)", prefix, idx[, "row"], idx[, "col"]))
}
