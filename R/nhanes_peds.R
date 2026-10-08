# Internal helpers for sample_nhanes_peds() and compare_nhanes_peds().

nhanes_required_cols <- c("CYCLE", "SEQN", "AGE", "SEX", "WT", "HT", "MEC_WT")

# Children eligible as donors (all `vars` measured), with sampling
# probabilities: MEC weight normalised within cycle x age x sex, times an
# equal share per cycle. Stops if any cycle x age x sex stratum is empty.
nhanes_donors <- function(data, ages, sex, vars, cycles) {
  keep <- data$AGE %in% ages & data$SEX %in% sex & data$CYCLE %in% cycles
  for (v in vars) {
    keep <- keep & !is.na(data[[v]])
  }
  dat <- data[keep, , drop = FALSE]

  grid <- expand.grid(
    CYCLE = cycles, AGE = ages, SEX = sex, stringsAsFactors = FALSE
  )
  empty <- !paste(grid$CYCLE, grid$AGE, grid$SEX) %in%
    paste(dat$CYCLE, dat$AGE, dat$SEX)
  if (any(empty)) {
    stop(
      "No eligible children for: ",
      paste(grid$CYCLE[empty], "age", grid$AGE[empty], grid$SEX[empty],
            collapse = "; "),
      call. = FALSE
    )
  }

  cycle_total <- stats::ave(dat$MEC_WT, dat$CYCLE, dat$AGE, dat$SEX, FUN = sum)
  dat$PROB <- dat$MEC_WT / cycle_total / length(cycles)
  dat
}

# Gaussian kernel on log(vars) for one age x sex stratum of donors.
# H = (bandwidth_factor * 0.9)^2 * n_eff^(-2 / (d + 4)) * S R S.
nhanes_kernel <- function(donors, vars, bandwidth_factor) {
  d <- length(vars)
  p <- donors$PROB / sum(donors$PROB)
  y <- log(as.matrix(donors[vars]))

  sigma <- vapply(seq_len(d), function(j) robust_scale(y[, j], p), numeric(1))
  R <- diag(d)
  rho <- NA_real_
  if (d == 2) {
    rho <- if (all(sigma > 0)) weighted_cor(y[, 1], y[, 2], p) else 0
    R[1, 2] <- R[2, 1] <- rho
  }

  neff <- 1 / sum(p^2)
  S <- diag(sigma, nrow = d)
  H <- (bandwidth_factor * 0.9)^2 * neff^(-2 / (d + 4)) * (S %*% R %*% S)

  list(
    y = y,
    p = p,
    H = H,
    neff = neff,
    rho = rho,
    cycle = donors$CYCLE,
    seqn = donors$SEQN
  )
}
