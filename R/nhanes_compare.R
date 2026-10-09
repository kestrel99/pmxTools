nhanes_summary <- function(x, w) {
  w <- w / sum(w)
  q <- weighted_quantile(x, w, c(0.05, 0.25, 0.5, 0.75, 0.95))
  c(
    N = length(x), Mean = sum(w * x), SD = weighted_sd(x, w),
    P05 = q[1], Q1 = q[2], Median = q[3], Q3 = q[4], P95 = q[5]
  )
}

#' Compare a simulated pediatric population with NHANES
#'
#' Summarises body weight, height and (when both were simulated) body mass
#' index by age and sex for a population from [sample_nhanes()] and for
#' the NHANES reference it was drawn from. The reference uses the same
#' releases and donor eligibility as the simulation, weighted with the
#' cycle-balanced MEC weights.
#'
#' @param sim Output of [sample_nhanes()], with its attributes intact.
#' @param data Reference data with the columns of [nhanes_peds]; normally the
#'   data used for the simulation.
#' @return A tibble with one row per `AGE`, `SEX` and `VARIABLE` (`WT`, `HT`
#'   and `BMI` in kg/m^2, as simulated). Columns `N`, `Mean`, `SD`, `P05`,
#'   `Q1`, `Median`, `Q3` and `P95` appear with suffix `_REF` (weighted
#'   NHANES reference) and `_SIM` (simulated), followed by the ratios
#'   `MEAN_RATIO`, `MEDIAN_RATIO`, `P05_RATIO` and `P95_RATIO` (simulated /
#'   reference). When both weight and height were simulated, `RHO_REF` and
#'   `RHO_SIM` give the correlation of log weight and log height.
#' @seealso [sample_nhanes()], [plot_nhanes()]
#' @examples
#' sim <- sample_nhanes(n = 200, ages = c(4, 12), seed = 1)
#' compare_nhanes(sim)
#' @export
compare_nhanes <- function(sim, data = pmxTools::nhanes_peds) {
  cycles <- attr(sim, "cycles")
  vars <- attr(sim, "vars")
  if (is.null(cycles) || is.null(vars)) {
    stop(
      "'sim' must be the output of sample_nhanes() ",
      "(its \"cycles\" and \"vars\" attributes are missing)",
      call. = FALSE
    )
  }
  missing_cols <- setdiff(c("AGE", "SEX", vars), names(sim))
  if (length(missing_cols) > 0) {
    stop("'sim' is missing column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  ages <- sort(unique(sim$AGE))
  sex <- unique(sim$SEX)
  ref <- nhanes_donors(data, ages, sex, vars, cycles)
  sim <- as.data.frame(sim)

  both <- length(vars) == 2
  out_vars <- vars
  if (both) {
    ref$BMI <- ref$WT / (ref$HT / 100)^2
    sim$BMI <- sim$WT / (sim$HT / 100)^2
    out_vars <- c(vars, "BMI")
  }

  rows <- list()
  for (age in ages) {
    for (sx in sex) {
      r <- ref[ref$AGE == age & ref$SEX == sx, , drop = FALSE]
      s <- sim[sim$AGE == age & sim$SEX == sx, , drop = FALSE]
      for (v in out_vars) {
        ref_stats <- nhanes_summary(r[[v]], r$PROB)
        sim_stats <- nhanes_summary(s[[v]], rep(1, nrow(s)))
        row <- data.frame(AGE = age, SEX = sx, VARIABLE = v,
                          stringsAsFactors = FALSE)
        row[paste0(names(ref_stats), "_REF")] <- as.list(ref_stats)
        row[paste0(names(sim_stats), "_SIM")] <- as.list(sim_stats)
        if (both) {
          row$RHO_REF <- weighted_cor(log(r$WT), log(r$HT), r$PROB)
          row$RHO_SIM <- stats::cor(log(s$WT), log(s$HT))
        }
        rows[[length(rows) + 1]] <- row
      }
    }
  }

  cmp <- do.call(rbind, rows)
  cmp$MEAN_RATIO <- cmp$Mean_SIM / cmp$Mean_REF
  cmp$MEDIAN_RATIO <- cmp$Median_SIM / cmp$Median_REF
  cmp$P05_RATIO <- cmp$P05_SIM / cmp$P05_REF
  cmp$P95_RATIO <- cmp$P95_SIM / cmp$P95_REF
  if (both) {
    rho <- c("RHO_REF", "RHO_SIM")
    cmp <- cmp[c(setdiff(names(cmp), rho), rho)]
  }
  tibble::as_tibble(cmp)
}

#' Plot simulated against NHANES percentiles
#'
#' Plots the 5th, 50th and 95th percentiles of each simulated measure against
#' age, for the simulated population and the NHANES reference, by sex.
#'
#' @param comparison Output of [compare_nhanes()].
#' @return A ggplot object.
#' @seealso [compare_nhanes()], [sample_nhanes()]
#' @examples
#' sim <- sample_nhanes(n = 200, ages = 2:17, seed = 1)
#' plot_nhanes(compare_nhanes(sim))
#' @export
plot_nhanes <- function(comparison) {
  percentiles <- c("P05", "Median", "P95")
  sources <- c(REF = "NHANES", SIM = "Simulated")
  long <- do.call(rbind, lapply(names(sources), function(src) {
    do.call(rbind, lapply(percentiles, function(pct) {
      data.frame(
        AGE = comparison$AGE,
        SEX = comparison$SEX,
        VARIABLE = comparison$VARIABLE,
        Percentile = pct,
        Source = sources[[src]],
        VALUE = comparison[[paste0(pct, "_", src)]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  long$Percentile <- factor(long$Percentile, levels = percentiles)
  long$VARIABLE <- factor(
    long$VARIABLE,
    levels = intersect(c("WT", "HT", "BMI"), long$VARIABLE)
  )

  ggplot(long, aes(x = .data$AGE, y = .data$VALUE,
                   colour = .data$Source, linetype = .data$Percentile)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1) +
    facet_grid(VARIABLE ~ SEX, scales = "free_y") +
    labs(
      title = "NHANES vs simulated percentiles",
      x = "Age (years)", y = NULL, colour = NULL
    ) +
    theme_bw()
}
