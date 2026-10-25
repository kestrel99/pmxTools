# Heterogeneous correlation matrix: Pearson (numeric-numeric), polyserial
# (numeric-ordinal) and polychoric (ordinal-ordinal). Polyserial/polychoric
# use the two-step ML estimator: thresholds from the marginal proportions,
# then ML for rho.

#' Heterogeneous correlation matrix
#'
#' Computes a correlation matrix for a mix of continuous and ordinal
#' variables, choosing the estimator for each pair: Pearson for two numeric
#' variables, polyserial for a numeric and an ordinal variable, and polychoric
#' for two ordinal variables. Useful for exploring covariate relationships
#' before covariate modelling or simulation.
#'
#' Polychoric and polyserial correlations estimate the correlation of
#' underlying normal variables from which the ordinal categories are assumed
#' to arise. Both use the two-step estimator: the thresholds between
#' categories are set from the marginal proportions, and the correlation is
#' then estimated by maximum likelihood (Olsson 1979; Olsson, Drasgow and
#' Dorans 1982).
#'
#' Numeric columns are treated as continuous and factor columns as ordinal,
#' using their level order. Unordered factors are accepted with a warning,
#' because their level order (often alphabetical) may not be meaningful;
#' convert them with `factor(x, levels = ..., ordered = TRUE)` first.
#' Each correlation uses the cases complete for that pair.
#'
#' @param data A data frame of numeric and/or factor columns.
#' @param plot If `TRUE`, also draw the matrix with [plot_cor_triangle()].
#' @param ... Passed to [plot_cor_triangle()] (e.g. `labels`, `text_size`).
#' @return A symmetric correlation matrix with the column names of `data`.
#'   Pairs that cannot be estimated (a factor with one level, or a constant
#'   numeric column) are `NA`. If `plot = TRUE`, a list with elements `cor`
#'   (the matrix) and `plot` (a ggplot).
#' @references
#' Olsson U (1979). Maximum likelihood estimation of the polychoric
#' correlation coefficient. *Psychometrika*, 44, 443-460.
#'
#' Olsson U, Drasgow F, Dorans NJ (1982). The polyserial correlation
#' coefficient. *Psychometrika*, 47, 337-347.
#' @seealso [polychoric()], [polyserial()], [plot_cor_triangle()]
#' @examples
#' set.seed(1)
#' n <- 300
#' wt <- rnorm(n, 70, 12)
#' age <- 30 + 0.4 * (wt - 70) + rnorm(n, 0, 8)
#' egfr <- cut(-0.03 * age + rnorm(n), c(-Inf, -2, -1, Inf),
#'             labels = c("severe", "moderate", "normal"),
#'             ordered_result = TRUE)
#' covs <- data.frame(WT = wt, AGE = age, RENAL = egfr)
#' het_cor(covs)
#'
#' res <- het_cor(covs, plot = TRUE, labels = c(WT = "Weight (kg)"))
#' res$plot
#' @export
het_cor <- function(data, plot = FALSE, ...) {
  types <- cor_types(data)
  warn_unordered_factors(data)
  vars <- names(data)
  R <- diag(length(vars))
  dimnames(R) <- list(vars, vars)

  for (j in seq_len(length(vars) - 1)) {
    for (i in (j + 1):length(vars)) {
      ok <- !is.na(data[[i]]) & !is.na(data[[j]])
      x <- data[[i]][ok]
      y <- data[[j]][ok]
      R[i, j] <- R[j, i] <- switch(types[i, j],
        Pearson = stats::cor(x, y),
        Polychoric = polychoric_rho(x, y),
        Polyserial = if (is.numeric(x)) polyserial_rho(x, y) else polyserial_rho(y, x)
      )
    }
  }

  if (!plot) return(R)
  list(cor = R, plot = plot_cor_triangle(R, types = types, ...))
}

# Correlation type for each pair of columns: "Pearson", "Polyserial" or
# "Polychoric", by how many of the two columns are ordinal (factor).
cor_types <- function(data) {
  if (!is.data.frame(data) || ncol(data) < 2) {
    stop("'data' must be a data frame with at least two columns", call. = FALSE)
  }
  is_num <- vapply(data, is.numeric, logical(1))
  is_ord <- vapply(data, is.factor, logical(1))
  bad <- names(data)[!(is_num | is_ord)]
  if (length(bad)) {
    stop("Columns must be numeric or factor: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  n_ord <- outer(is_ord, is_ord, "+")
  matrix(c("Pearson", "Polyserial", "Polychoric")[n_ord + 1], ncol(data),
         dimnames = list(names(data), names(data)))
}

warn_unordered_factors <- function(data) {
  unordered <- names(data)[vapply(data, function(v) {
    is.factor(v) && !is.ordered(v)
  }, logical(1))]
  if (length(unordered)) {
    warning(
      "Treating unordered factor(s) as ordinal in level order: ",
      paste(unordered, collapse = ", "),
      call. = FALSE
    )
  }
}

warn_unordered <- function(x, arg) {
  if (!is.ordered(x)) {
    warning("Treating unordered factor '", arg, "' as ordinal in level order",
            call. = FALSE)
  }
}

# Normal thresholds implied by the cumulative marginal proportions.
cut_points <- function(counts) {
  stats::qnorm(cumsum(counts) / sum(counts))[-length(counts)]
}

# Bivariate standard normal CDF P(X <= h, Y <= k) with correlation rho,
# via P = integral_{-Inf}^{h} phi(u) * Phi((k - rho * u) / sqrt(1 - rho^2)) du.
pbvnorm <- function(h, k, rho) {
  if (h == -Inf || k == -Inf) return(0)
  if (h == Inf) return(stats::pnorm(k))
  if (k == Inf) return(stats::pnorm(h))
  s <- sqrt(1 - rho^2)
  f <- function(u) stats::dnorm(u) * stats::pnorm((k - rho * u) / s)
  stats::integrate(f, -Inf, h, rel.tol = 1e-10)$value
}

polychoric_rho <- function(x, y, max_cor = 0.9999) {
  tab <- table(droplevels(x), droplevels(y))
  if (any(dim(tab) < 2)) return(NA_real_)
  a <- c(-Inf, cut_points(rowSums(tab)), Inf)
  b <- c(-Inf, cut_points(colSums(tab)), Inf)

  nll <- function(rho) {
    cdf <- outer(a, b, Vectorize(function(h, k) pbvnorm(h, k, rho)))
    probs <- diff(t(diff(cdf)))  # cell probabilities, transposed
    -sum(t(tab) * log(pmax(probs, .Machine$double.eps)))
  }
  stats::optimize(nll, c(-max_cor, max_cor))$minimum
}

polyserial_rho <- function(x, y, max_cor = 0.9999) {
  y <- droplevels(y)
  if (nlevels(y) < 2 || length(x) < 2 || stats::sd(x) == 0) return(NA_real_)
  z <- (x - mean(x)) / stats::sd(x)
  tau <- c(-Inf, cut_points(table(y)), Inf)
  k <- as.integer(y)

  # Only P(y | z) depends on rho; the density of z drops out of the ML.
  nll <- function(rho) {
    s <- sqrt(1 - rho^2)
    p <- stats::pnorm((tau[k + 1] - rho * z) / s) -
      stats::pnorm((tau[k] - rho * z) / s)
    -sum(log(pmax(p, .Machine$double.eps)))
  }
  stats::optimize(nll, c(-max_cor, max_cor))$minimum
}

check_same_length <- function(x, y) {
  if (length(x) != length(y)) {
    stop("'x' and 'y' must be the same length", call. = FALSE)
  }
}

check_max_cor <- function(max_cor) {
  if (!is.numeric(max_cor) || length(max_cor) != 1 || is.na(max_cor) ||
        max_cor <= 0 || max_cor >= 1) {
    stop("'max_cor' must be a single number between 0 and 1 (exclusive)",
         call. = FALSE)
  }
}

#' Polychoric correlation
#'
#' Estimates the correlation between two ordinal variables, assuming each
#' arises from cutting an underlying normal variable into categories. Uses
#' the two-step estimator: thresholds from the marginal proportions, then
#' maximum likelihood for the correlation (Olsson 1979).
#'
#' @param x,y Factors of the same length. Their level order is used as the
#'   ordinal scale; unordered factors give a warning.
#' @param max_cor Largest absolute correlation searched, between 0 and 1.
#' @return The estimated correlation, using cases complete for both
#'   variables, or `NA` if either has fewer than two observed levels.
#' @references
#' Olsson U (1979). Maximum likelihood estimation of the polychoric
#' correlation coefficient. *Psychometrika*, 44, 443-460.
#' @seealso [het_cor()], [polyserial()]
#' @examples
#' set.seed(1)
#' z1 <- rnorm(500)
#' z2 <- 0.6 * z1 + sqrt(1 - 0.6^2) * rnorm(500)
#' x <- cut(z1, c(-Inf, -0.5, 0.5, Inf), ordered_result = TRUE)
#' y <- cut(z2, c(-Inf, 0, 1, Inf), ordered_result = TRUE)
#' polychoric(x, y)
#' @export
polychoric <- function(x, y, max_cor = 0.9999) {
  if (!is.factor(x) || !is.factor(y)) {
    stop("'x' and 'y' must be factors", call. = FALSE)
  }
  check_same_length(x, y)
  check_max_cor(max_cor)
  warn_unordered(x, "x")
  warn_unordered(y, "y")
  ok <- !is.na(x) & !is.na(y)
  polychoric_rho(x[ok], y[ok], max_cor)
}

#' Polyserial correlation
#'
#' Estimates the correlation between a continuous variable and an ordinal
#' variable, assuming the ordinal variable arises from cutting an underlying
#' normal variable into categories. Uses the two-step estimator: thresholds
#' from the marginal proportions, then maximum likelihood for the
#' correlation (Olsson, Drasgow and Dorans 1982).
#'
#' @param x A numeric vector.
#' @param y A factor of the same length as `x`. Its level order is used as the
#'   ordinal scale; an unordered factor gives a warning.
#' @param max_cor Largest absolute correlation searched, between 0 and 1.
#' @return The estimated correlation, using cases complete for both
#'   variables, or `NA` if `y` has fewer than two observed levels or `x` is
#'   constant.
#' @references
#' Olsson U, Drasgow F, Dorans NJ (1982). The polyserial correlation
#' coefficient. *Psychometrika*, 47, 337-347.
#' @seealso [het_cor()], [polychoric()]
#' @examples
#' set.seed(1)
#' x <- rnorm(500)
#' latent <- 0.6 * x + sqrt(1 - 0.6^2) * rnorm(500)
#' y <- cut(latent, c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
#' polyserial(x, y)
#' @export
polyserial <- function(x, y, max_cor = 0.9999) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric", call. = FALSE)
  }
  if (!is.factor(y)) {
    stop("'y' must be a factor", call. = FALSE)
  }
  check_same_length(x, y)
  check_max_cor(max_cor)
  warn_unordered(y, "y")
  ok <- !is.na(x) & !is.na(y)
  polyserial_rho(x[ok], y[ok], max_cor)
}
