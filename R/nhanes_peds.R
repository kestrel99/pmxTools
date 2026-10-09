# Internal helpers for sample_nhanes() and compare_nhanes().

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

# Draw n children from one stratum's kernel.
nhanes_sample_kernel <- function(kernel, n, vars) {
  d <- length(vars)
  idx <- sample.int(length(kernel$p), size = n, replace = TRUE, prob = kernel$p)
  noise <- matrix(MASS::mvrnorm(n, mu = rep(0, d), Sigma = kernel$H), nrow = n)
  values <- exp(kernel$y[idx, , drop = FALSE] + noise)
  colnames(values) <- vars
  out <- as.data.frame(values)
  out$SOURCE_CYCLE <- kernel$cycle[idx]
  out$SOURCE_SEQN <- kernel$seqn[idx]
  out
}

nhanes_check_values <- function(arg, values, allowed) {
  bad <- setdiff(values, allowed)
  if (length(bad) > 0) {
    stop(
      "'", arg, "' value(s) not in the reference data: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
}

nhanes_check_vars <- function(vars) {
  if (!is.character(vars) || length(vars) == 0 ||
      !all(vars %in% c("WT", "HT"))) {
    stop("'vars' must be \"WT\", \"HT\" or both", call. = FALSE)
  }
  intersect(c("WT", "HT"), vars)
}

#' Simulate pediatric body weight and height from NHANES
#'
#' Simulates a virtual pediatric population with realistic body weight and/or
#' standing height by smoothed resampling of children in the National Health
#' and Nutrition Examination Survey (NHANES), see [nhanes_peds].
#'
#' Each age (whole years) x sex stratum is simulated separately:
#'
#' 1. Donors are the stratum's children with every requested measure in
#'    `vars`. A child missing height is therefore still used when only weight
#'    is requested.
#' 2. Each donor's sampling probability is its MEC exam weight, normalised
#'    within its NHANES release, times an equal share per release, so that no
#'    release dominates because of a larger total weight. The equal share is a
#'    modelling choice, not an official pooled NHANES weight.
#' 3. `n` donors are drawn with these probabilities and Gaussian noise is
#'    added on the log scale, so values stay positive and continuous. The
#'    noise covariance is
#'    \deqn{H = (b \cdot 0.9)^2 \, n_{eff}^{-2/(d+4)} \, S R S}
#'    where `b` is `bandwidth_factor`, `d` is the number of `vars`,
#'    \eqn{n_{eff}} is Kish's effective sample size, `S` holds each log
#'    measure's robust scale \eqn{\min(SD, IQR/1.349)} and `R` is the weighted
#'    correlation matrix of the log measures. With one measure this is the
#'    robust Silverman rule \eqn{h = 0.9 \, \sigma \, n_{eff}^{-1/5}}; with
#'    both, the correlated noise preserves the weight-height relationship.
#'
#' Use [compare_nhanes()] and [plot_nhanes()] to check the
#' simulated population against the reference.
#'
#' @param n Number of children to simulate per age x sex stratum.
#' @param ages Ages in whole years (2 to 17 for [nhanes_peds]).
#' @param sex `"Male"`, `"Female"` or both.
#' @param vars Measures to simulate: `"WT"` (body weight, kg), `"HT"`
#'   (standing height, cm) or both.
#' @param cycles NHANES releases to use (values of `data$CYCLE`); `NULL` uses
#'   all of them.
#' @param bandwidth_factor Multiplier for the kernel bandwidth. `0` gives plain
#'   weighted resampling of real children; values above 1 smooth more.
#' @param seed Optional random seed. The previous random number generator
#'   state is restored afterwards.
#' @param data Reference data with the columns of [nhanes_peds].
#' @return A tibble with one row per simulated child and columns `ID`, `AGE`,
#'   `SEXN` (1 = male, 2 = female), `SEX`, the requested `vars`,
#'   `SOURCE_CYCLE` and `SOURCE_SEQN` (the NHANES release and sequence number
#'   of the donor child). Attributes: `"kernels"`, a tibble of per-stratum
#'   diagnostics (`N_NHANES` donors, `N_EFF` effective sample size, `H_WT`
#'   and/or `H_HT` bandwidth SDs on the log scale, and `RHO`, the weighted
#'   log weight-height correlation, when both measures are simulated);
#'   `"cycles"` and `"vars"`, as used.
#' @seealso [nhanes_peds], [compare_nhanes()], [plot_nhanes()]
#' @examples
#' sim <- sample_nhanes(n = 100, ages = c(2, 8, 14), seed = 20261008)
#' head(sim)
#' attr(sim, "kernels")
#'
#' # Weight only, using every child with a measured weight
#' wt <- sample_nhanes(n = 100, ages = 10, vars = "WT", seed = 1)
#' summary(wt$WT)
#' @export
sample_nhanes <- function(n = 500,
                               ages = 2:17,
                               sex = c("Male", "Female"),
                               vars = c("WT", "HT"),
                               cycles = NULL,
                               bandwidth_factor = 1,
                               seed = NULL,
                               data = pmxTools::nhanes_peds) {
  missing_cols <- setdiff(nhanes_required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "'data' is missing column(s): ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 || n != round(n)) {
    stop("'n' must be a single positive whole number", call. = FALSE)
  }
  if (!is.numeric(bandwidth_factor) || length(bandwidth_factor) != 1 ||
      is.na(bandwidth_factor) || bandwidth_factor < 0) {
    stop("'bandwidth_factor' must be a single non-negative number", call. = FALSE)
  }
  vars <- nhanes_check_vars(vars)
  if (is.null(cycles)) {
    cycles <- sort(unique(data$CYCLE))
  }
  nhanes_check_values("ages", ages, data$AGE)
  nhanes_check_values("sex", sex, data$SEX)
  nhanes_check_values("cycles", cycles, data$CYCLE)
  ages <- sort(unique(ages))
  sex <- unique(sex)

  if (!is.null(seed)) {
    genv <- globalenv()
    old_seed <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      get(".Random.seed", envir = genv, inherits = FALSE)
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = genv)
      } else {
        assign(".Random.seed", old_seed, envir = genv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  donors <- nhanes_donors(data, ages, sex, vars, cycles)
  strata <- expand.grid(SEX = sex, AGE = ages, stringsAsFactors = FALSE)

  sims <- vector("list", nrow(strata))
  kernels <- vector("list", nrow(strata))
  for (i in seq_len(nrow(strata))) {
    age <- strata$AGE[i]
    sx <- strata$SEX[i]
    stratum <- donors[donors$AGE == age & donors$SEX == sx, , drop = FALSE]
    kernel <- nhanes_kernel(stratum, vars, bandwidth_factor)

    sims[[i]] <- cbind(
      data.frame(AGE = age, SEX = sx, stringsAsFactors = FALSE),
      nhanes_sample_kernel(kernel, n, vars)
    )

    diagnostics <- data.frame(
      AGE = age, SEX = sx, N_NHANES = nrow(stratum), N_EFF = kernel$neff,
      stringsAsFactors = FALSE
    )
    for (j in seq_along(vars)) {
      diagnostics[[paste0("H_", vars[j])]] <- sqrt(kernel$H[j, j])
    }
    if (length(vars) == 2) {
      diagnostics$RHO <- kernel$rho
    }
    kernels[[i]] <- diagnostics
  }

  sim <- do.call(rbind, sims)
  sim$ID <- seq_len(nrow(sim))
  sim$SEXN <- ifelse(sim$SEX == "Male", 1L, 2L)
  sim <- tibble::as_tibble(
    sim[c("ID", "AGE", "SEXN", "SEX", vars, "SOURCE_CYCLE", "SOURCE_SEQN")]
  )

  attr(sim, "kernels") <- tibble::as_tibble(do.call(rbind, kernels))
  attr(sim, "cycles") <- cycles
  attr(sim, "vars") <- vars
  sim
}
