# Internal weighted statistics used by sample_nhanes_peds() and
# compare_nhanes_peds(). Weights need not sum to 1.

weighted_quantile <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w)
  stats::approx(
    x = c(0, cumsum(w)),
    y = c(x[1], x),
    xout = p,
    ties = "ordered",
    rule = 2
  )$y
}

weighted_sd <- function(x, w) {
  w <- w / sum(w)
  mu <- sum(w * x)
  sqrt(sum(w * (x - mu)^2))
}

weighted_cor <- function(x, y, w) {
  w <- w / sum(w)
  dx <- x - sum(w * x)
  dy <- y - sum(w * y)
  sum(w * dx * dy) / sqrt(sum(w * dx^2) * sum(w * dy^2))
}

# Robust Silverman-type scale: min(SD, IQR / 1.349), falling back to the SD
# and then to 0 for degenerate data.
robust_scale <- function(x, w) {
  s <- weighted_sd(x, w)
  iqr <- diff(weighted_quantile(x, w, c(0.25, 0.75))) / 1.349
  scale <- min(s, iqr)
  if (!is.finite(scale) || scale <= 0) {
    scale <- s
  }
  if (!is.finite(scale) || scale <= 0) {
    scale <- 0
  }
  scale
}
