# NHANES Pediatric Weight/Height Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `sample_nhanes_peds()`, `compare_nhanes_peds()`, `plot_nhanes_peds()` and the bundled `nhanes_peds` dataset to pmxTools, turning the stand-alone NHANES script into package functionality that simulates pediatric weight and/or height.

**Architecture:** A data-raw script builds a small bundled dataset from four NHANES cycles (ages 2-17, weight and height, `NA` kept). Internal weighted-statistics helpers feed a kernel builder that, per age/sex stratum, resamples cycle-balanced MEC-weighted donors and adds 1-D or 2-D Gaussian noise on the log scale (Silverman-type bandwidth). Validation is separate: a comparison table against the weighted reference and a ggplot.

**Tech Stack:** R (>= 4.0), base R + `MASS::mvrnorm`, `tibble`, `ggplot2` (all already imported); `haven` and `usethis` only in `data-raw/`; testthat 3e.

**Spec:** `dev/specs/2026-10-08-nhanes-peds-design.md`. Branch: `nhanes-peds`.

**Conventions for every task**

- Run commands from the package root `C:/Users/justin/Documents/GitHub/pmxTools` (Git Bash).
- Run a single test file with:
  `Rscript -e "devtools::test(filter = '<name>', reporter = 'summary')"`
- On Windows, running tests deletes the vdiffr snapshots under `tests/testthat/_snaps/plot/` and rewrites line endings in `_snaps/*.md`. After every test run, before committing: `git checkout -- tests/testthat/_snaps`.
- Internal (non-exported) functions are visible to tests because `devtools::test()` loads the namespace.
- Error messages use `stop(..., call. = FALSE)`.

---

## File structure

| File | Status | Responsibility |
|------|--------|----------------|
| `data-raw/nhanes_peds.R` | create | Download NHANES cycles, build and save `nhanes_peds` |
| `data/nhanes_peds.rda` | generated | Bundled reference data |
| `R/data.R` | create | Roxygen docs for `nhanes_peds` |
| `R/weighted_stats.R` | create | Internal weighted quantile, SD, correlation, robust scale |
| `R/nhanes_peds.R` | create | `sample_nhanes_peds()` + internal donor selection, kernel, sampling |
| `R/nhanes_compare.R` | create | `compare_nhanes_peds()`, `plot_nhanes_peds()` |
| `tests/testthat/test-weighted_stats.R` | create | Helper tests |
| `tests/testthat/test-nhanes_peds.R` | create | Data + simulator tests |
| `tests/testthat/test-nhanes_compare.R` | create | Comparison + plot tests |
| `DESCRIPTION` | modify | `LazyData: true` (added by `usethis::use_data()`) |
| `NEWS.md` | modify | Development-version bullet |

---

### Task 1: Bundled dataset `nhanes_peds`

**Files:**
- Create: `data-raw/nhanes_peds.R`
- Create (generated): `data/nhanes_peds.rda`
- Create: `R/data.R`
- Create: `tests/testthat/test-nhanes_peds.R`
- Modify: `DESCRIPTION` (via `usethis::use_data()`)

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-nhanes_peds.R`:

```r
test_that("nhanes_peds has the documented structure", {
  expect_s3_class(nhanes_peds, "data.frame")
  expect_named(
    nhanes_peds,
    c("CYCLE", "SEQN", "AGE_MONTHS", "AGE", "SEX", "WT", "HT", "MEC_WT")
  )
  expect_setequal(unique(nhanes_peds$CYCLE), c("2013-14", "2015-16", "2017-18", "2021-23"))
  expect_setequal(unique(nhanes_peds$AGE), 2:17)
  expect_setequal(unique(nhanes_peds$SEX), c("Male", "Female"))
  expect_false(anyNA(nhanes_peds$MEC_WT))
  expect_true(all(nhanes_peds$MEC_WT > 0))
})

test_that("nhanes_peds keeps children with only one of weight or height", {
  has_wt <- !is.na(nhanes_peds$WT)
  has_ht <- !is.na(nhanes_peds$HT)
  expect_true(all(has_wt | has_ht))
  expect_true(any(has_wt & !has_ht))
  expect_true(all(nhanes_peds$WT[has_wt] > 0))
  expect_true(all(nhanes_peds$HT[has_ht] > 0))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: errors `object 'nhanes_peds' not found`.

- [ ] **Step 3: Write the data-raw script**

Create `data-raw/nhanes_peds.R`:

```r
# Build data/nhanes_peds.rda from the public NHANES files.
#
# Run from the package root:  Rscript data-raw/nhanes_peds.R
# Needs internet access and the haven and usethis packages.
#
# Children aged 2-17 whole years at the examination, from four NHANES
# releases. Children are kept when at least one of body weight or standing
# height was measured; a missing measure is kept as NA.

cycles <- data.frame(
  CYCLE = c("2013-14", "2015-16", "2017-18", "2021-23"),
  YEAR = c(2013, 2015, 2017, 2021),
  SUFFIX = c("H", "I", "J", "L")
)

read_nhanes_cycle <- function(cycle, year, suffix) {
  url <- paste0(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", year, "/DataFiles/"
  )
  demo <- haven::read_xpt(paste0(url, "DEMO_", suffix, ".XPT"))
  bmx <- haven::read_xpt(paste0(url, "BMX_", suffix, ".XPT"))

  demo_vars <- c("SEQN", "RIAGENDR", "RIDEXAGM", "WTMEC2YR")
  bmx_vars <- c("SEQN", "BMXWT", "BMXHT")
  missing_vars <- c(setdiff(demo_vars, names(demo)), setdiff(bmx_vars, names(bmx)))
  if (length(missing_vars) > 0) {
    stop("Cycle ", cycle, " lacks: ", paste(missing_vars, collapse = ", "))
  }

  d <- merge(demo[demo_vars], bmx[bmx_vars], by = "SEQN")
  positive <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x) | x <= 0] <- NA_real_
    x
  }

  out <- data.frame(
    CYCLE = cycle,
    SEQN = as.numeric(d$SEQN),
    AGE_MONTHS = as.numeric(d$RIDEXAGM),
    AGE = as.integer(floor(as.numeric(d$RIDEXAGM) / 12)),
    SEX = c("Male", "Female")[match(as.integer(d$RIAGENDR), 1:2)],
    WT = positive(d$BMXWT),
    HT = positive(d$BMXHT),
    MEC_WT = positive(d$WTMEC2YR),
    stringsAsFactors = FALSE
  )

  out[
    !is.na(out$AGE) & out$AGE %in% 2:17 &
      !is.na(out$SEX) &
      !is.na(out$MEC_WT) &
      (!is.na(out$WT) | !is.na(out$HT)),
  ]
}

nhanes_peds <- do.call(
  rbind,
  Map(read_nhanes_cycle, cycles$CYCLE, cycles$YEAR, cycles$SUFFIX)
)
nhanes_peds <- nhanes_peds[order(nhanes_peds$CYCLE, nhanes_peds$SEQN), ]
rownames(nhanes_peds) <- NULL

message("Rows: ", nrow(nhanes_peds))
message("Missing WT only: ", sum(is.na(nhanes_peds$WT)))
message("Missing HT only: ", sum(is.na(nhanes_peds$HT)))

usethis::use_data(nhanes_peds, overwrite = TRUE, compress = "xz")
```

- [ ] **Step 4: Build the data**

Run: `Rscript data-raw/nhanes_peds.R`
Expected: messages `Rows: ...` (about 11,300), `Missing WT only: ...`, `Missing HT only: ...` (both > 0), then usethis output saving `data/nhanes_peds.rda` and setting `LazyData: true` in `DESCRIPTION`. Record the three numbers in the commit message.

Check: `git diff DESCRIPTION` shows `LazyData: true` (and possibly `Depends: R (>= 3.5)` — if usethis lowered the existing `R (>= 4.0.0)`, restore `R (>= 4.0.0)`).

- [ ] **Step 5: Document the dataset**

Create `R/data.R`:

```r
#' NHANES pediatric body weight and height reference data
#'
#' Body weight and standing height of children aged 2 to 17 years at the
#' examination, from four releases of the US National Health and Nutrition
#' Examination Survey (NHANES): 2013-2014, 2015-2016, 2017-2018 and
#' August 2021-August 2023. Used as the reference population by
#' [sample_nhanes_peds()].
#'
#' Children are included when they have a positive two-year mobile
#' examination center (MEC) exam weight and at least one of body weight or
#' standing height. A child with only one of the two measures is kept, with
#' the other set to `NA`; [sample_nhanes_peds()] uses such children only when
#' the missing measure is not requested.
#'
#' The MEC weights are those of each individual release. They are not
#' combined into an official pooled weight; [sample_nhanes_peds()] gives each
#' release an equal share, which is a modelling choice.
#'
#' NHANES data are produced by the US National Center for Health Statistics
#' and are in the public domain. The data were built with
#' `data-raw/nhanes_peds.R` in the package source.
#'
#' @format A data frame with one row per child and 8 columns:
#' \describe{
#'   \item{CYCLE}{NHANES release, e.g. `"2017-18"`.}
#'   \item{SEQN}{NHANES respondent sequence number (unique within a release).}
#'   \item{AGE_MONTHS}{Age in months at the examination (`RIDEXAGM`).}
#'   \item{AGE}{Age in whole years, `floor(AGE_MONTHS / 12)`.}
#'   \item{SEX}{`"Male"` or `"Female"` (`RIAGENDR`).}
#'   \item{WT}{Body weight in kg (`BMXWT`); `NA` if not measured.}
#'   \item{HT}{Standing height in cm (`BMXHT`); `NA` if not measured.}
#'   \item{MEC_WT}{Two-year MEC exam weight (`WTMEC2YR`).}
#' }
#' @source National Center for Health Statistics, NHANES public data files
#'   `DEMO_H`/`BMX_H`, `DEMO_I`/`BMX_I`, `DEMO_J`/`BMX_J` and
#'   `DEMO_L`/`BMX_L`, \url{https://wwwn.cdc.gov/nchs/nhanes/}.
#' @seealso [sample_nhanes_peds()]
"nhanes_peds"
```

Run: `Rscript -e "devtools::document()"`
Expected: `Writing 'nhanes_peds.Rd'`. (A warning that `sample_nhanes_peds` is an unresolved link is expected until Task 4.)

- [ ] **Step 6: Run test to verify it passes**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: both tests pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 7: Commit**

Use the three numbers printed in Step 4 in the message body, e.g.
`Rows: 11300; missing WT: 12; missing HT: 150.`

```bash
git add data-raw/nhanes_peds.R data/nhanes_peds.rda R/data.R man/nhanes_peds.Rd DESCRIPTION tests/testthat/test-nhanes_peds.R
git commit -m "Add nhanes_peds reference data (NHANES 2013-2023, ages 2-17)" -m "Rows: R; missing WT: W; missing HT: H."   # substitute R, W, H from Step 4
```

---

### Task 2: Weighted statistics helpers

**Files:**
- Create: `R/weighted_stats.R`
- Create: `tests/testthat/test-weighted_stats.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-weighted_stats.R`:

```r
test_that("weighted_quantile matches the unweighted median for equal weights", {
  x <- c(5, 1, 3, 2, 4)
  expect_equal(weighted_quantile(x, rep(1, 5), 0.5), 3)
})

test_that("weighted_quantile follows the weights", {
  x <- c(1, 2, 3)
  expect_equal(weighted_quantile(x, c(0, 0, 1), 0.5), 3)
  expect_equal(weighted_quantile(x, c(1, 0, 0), c(0, 1)), c(1, 1))
})

test_that("weighted_sd is the population SD for equal weights", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  expect_equal(weighted_sd(x, rep(1, 8)), 2)
})

test_that("weighted_cor matches cor() for equal weights", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 1, 4, 3, 5)
  expect_equal(weighted_cor(x, y, rep(1, 5)), cor(x, y))
})

test_that("robust_scale is min(SD, IQR/1.349) and falls back to 0", {
  x <- c(1, 2, 3, 4, 100)
  w <- rep(1, 5)
  iqr <- diff(weighted_quantile(x, w, c(0.25, 0.75))) / 1.349
  expect_equal(robust_scale(x, w), min(weighted_sd(x, w), iqr))
  expect_equal(robust_scale(c(3, 3, 3), rep(1, 3)), 0)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'weighted_stats', reporter = 'summary')"`
Expected: errors `could not find function "weighted_quantile"` (and the others).

- [ ] **Step 3: Implement the helpers**

Create `R/weighted_stats.R`:

```r
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'weighted_stats', reporter = 'summary')"`
Expected: all 5 tests pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 5: Commit**

```bash
git add R/weighted_stats.R tests/testthat/test-weighted_stats.R
git commit -m "Add internal weighted statistics helpers"
```

---

### Task 3: Donor selection and kernel construction (internal)

**Files:**
- Create: `R/nhanes_peds.R`
- Modify: `tests/testthat/test-nhanes_peds.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-nhanes_peds.R`:

```r
test_that("nhanes_donors gives each cycle an equal share within a stratum", {
  donors <- nhanes_donors(
    nhanes_peds, ages = 8, sex = "Female", vars = "WT",
    cycles = unique(nhanes_peds$CYCLE)
  )
  shares <- tapply(donors$PROB, donors$CYCLE, sum)
  expect_equal(unname(shares), rep(0.25, 4))
})

test_that("nhanes_donors requires all requested vars", {
  cyc <- unique(nhanes_peds$CYCLE)
  wt_only <- nhanes_donors(nhanes_peds, 2:17, c("Male", "Female"), "WT", cyc)
  both <- nhanes_donors(nhanes_peds, 2:17, c("Male", "Female"), c("WT", "HT"), cyc)
  expect_false(anyNA(wt_only$WT))
  expect_true(anyNA(wt_only$HT))
  expect_false(anyNA(both$WT) || anyNA(both$HT))
  expect_gt(nrow(wt_only), nrow(both))
})

test_that("nhanes_donors stops on empty strata", {
  dat <- nhanes_peds[!(nhanes_peds$CYCLE == "2013-14" &
                         nhanes_peds$AGE == 5 &
                         nhanes_peds$SEX == "Male"), ]
  expect_error(
    nhanes_donors(dat, 5, "Male", "WT", unique(dat$CYCLE)),
    "No eligible children for: 2013-14 age 5 Male"
  )
})

test_that("1-D kernel bandwidth matches the original script", {
  donors <- nhanes_donors(
    nhanes_peds, ages = 10, sex = "Male", vars = "WT",
    cycles = unique(nhanes_peds$CYCLE)
  )
  k <- nhanes_kernel(donors, vars = "WT", bandwidth_factor = 1)

  # Script: h = 0.9 * min(sd, iqr / 1.349) * neff^(-1/5) on log(WT)
  p <- donors$PROB / sum(donors$PROB)
  x <- log(donors$WT)
  scale <- min(
    weighted_sd(x, p),
    diff(weighted_quantile(x, p, c(0.25, 0.75))) / 1.349
  )
  h <- 0.9 * scale * (1 / sum(p^2))^(-1 / 5)
  expect_equal(sqrt(k$H[1, 1]), h)
})

test_that("2-D kernel uses the d = 2 exponent and the weighted correlation", {
  donors <- nhanes_donors(
    nhanes_peds, ages = 10, sex = "Male", vars = c("WT", "HT"),
    cycles = unique(nhanes_peds$CYCLE)
  )
  k <- nhanes_kernel(donors, vars = c("WT", "HT"), bandwidth_factor = 1)
  p <- donors$PROB / sum(donors$PROB)
  neff <- 1 / sum(p^2)
  s_wt <- robust_scale(log(donors$WT), p)
  s_ht <- robust_scale(log(donors$HT), p)
  r <- weighted_cor(log(donors$WT), log(donors$HT), p)
  f <- 0.81 * neff^(-1 / 3)
  expect_equal(k$H, f * matrix(c(s_wt^2, r * s_wt * s_ht, r * s_wt * s_ht, s_ht^2), 2))
  expect_equal(k$rho, r)
  expect_gt(r, 0.5)
})

test_that("bandwidth_factor = 0 gives a zero bandwidth", {
  donors <- nhanes_donors(
    nhanes_peds, ages = 10, sex = "Male", vars = c("WT", "HT"),
    cycles = unique(nhanes_peds$CYCLE)
  )
  k <- nhanes_kernel(donors, vars = c("WT", "HT"), bandwidth_factor = 0)
  expect_equal(k$H, matrix(0, 2, 2))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: the two Task 1 tests pass; the new ones error `could not find function "nhanes_donors"` / `"nhanes_kernel"`.

- [ ] **Step 3: Implement donor selection and the kernel**

Create `R/nhanes_peds.R`:

```r
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: all tests pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 5: Commit**

```bash
git add R/nhanes_peds.R tests/testthat/test-nhanes_peds.R
git commit -m "Add NHANES donor selection and kernel construction"
```

---

### Task 4: `sample_nhanes_peds()` — interface, validation, seed

**Files:**
- Modify: `R/nhanes_peds.R` (append)
- Modify: `tests/testthat/test-nhanes_peds.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-nhanes_peds.R`:

```r
test_that("sample_nhanes_peds returns one block per age and sex", {
  sim <- sample_nhanes_peds(n = 7, ages = c(3, 12), seed = 1)
  expect_s3_class(sim, "tbl_df")
  expect_named(
    sim,
    c("ID", "AGE", "SEXN", "SEX", "WT", "HT", "SOURCE_CYCLE", "SOURCE_SEQN")
  )
  expect_equal(nrow(sim), 7 * 2 * 2)
  expect_equal(sim$ID, seq_len(28))
  expect_equal(as.vector(table(sim$AGE, sim$SEX)), rep(7L, 4))
  expect_equal(sim$SEXN, ifelse(sim$SEX == "Male", 1L, 2L))
  expect_true(all(is.finite(sim$WT) & sim$WT > 0))
  expect_true(all(is.finite(sim$HT) & sim$HT > 0))
})

test_that("sample_nhanes_peds stores kernels, cycles and vars", {
  sim <- sample_nhanes_peds(n = 5, ages = 4, sex = "Female", seed = 1)
  kern <- attr(sim, "kernels")
  expect_named(kern, c("AGE", "SEX", "N_NHANES", "N_EFF", "H_WT", "H_HT", "RHO"))
  expect_equal(nrow(kern), 1)
  expect_equal(attr(sim, "cycles"), sort(unique(nhanes_peds$CYCLE)))
  expect_equal(attr(sim, "vars"), c("WT", "HT"))
})

test_that("vars = 'WT' returns weight only", {
  sim <- sample_nhanes_peds(n = 5, ages = 4, vars = "WT", seed = 1)
  expect_named(
    sim, c("ID", "AGE", "SEXN", "SEX", "WT", "SOURCE_CYCLE", "SOURCE_SEQN")
  )
  expect_named(attr(sim, "kernels"), c("AGE", "SEX", "N_NHANES", "N_EFF", "H_WT"))
})

test_that("vars are returned in a fixed order", {
  sim <- sample_nhanes_peds(n = 2, ages = 4, vars = c("HT", "WT"), seed = 1)
  expect_equal(attr(sim, "vars"), c("WT", "HT"))
})

test_that("seed makes results reproducible without touching the global RNG", {
  set.seed(42)
  before <- get(".Random.seed", envir = globalenv())
  a <- sample_nhanes_peds(n = 5, ages = 6, seed = 123)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  b <- sample_nhanes_peds(n = 5, ages = 6, seed = 123)
  expect_identical(a, b)
})

test_that("cycles restricts the donors", {
  sim <- sample_nhanes_peds(n = 50, ages = 9, cycles = "2017-18", seed = 1)
  expect_true(all(sim$SOURCE_CYCLE == "2017-18"))
  expect_equal(attr(sim, "cycles"), "2017-18")
})

test_that("sample_nhanes_peds validates its arguments", {
  expect_error(sample_nhanes_peds(n = 0), "'n' must be a single positive whole number")
  expect_error(sample_nhanes_peds(n = 2.5), "'n' must be a single positive whole number")
  expect_error(sample_nhanes_peds(bandwidth_factor = -1), "'bandwidth_factor'")
  expect_error(sample_nhanes_peds(vars = "BMI"), "'vars' must be")
  expect_error(sample_nhanes_peds(vars = character(0)), "'vars' must be")
  expect_error(sample_nhanes_peds(ages = c(1, 18)), "'ages' value\\(s\\) not in the reference data: 1, 18")
  expect_error(sample_nhanes_peds(sex = "M"), "'sex' value\\(s\\) not in the reference data: M")
  expect_error(sample_nhanes_peds(cycles = "1999-00"), "'cycles' value\\(s\\) not in the reference data: 1999-00")
  expect_error(
    sample_nhanes_peds(data = nhanes_peds[c("AGE", "SEX")]),
    "'data' is missing column\\(s\\): CYCLE, SEQN, WT, HT, MEC_WT"
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: new tests error `could not find function "sample_nhanes_peds"`.

- [ ] **Step 3: Implement the sampler and the exported function**

Append to `R/nhanes_peds.R`:

```r
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
#' Use [compare_nhanes_peds()] and [plot_nhanes_peds()] to check the
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
#' @seealso [nhanes_peds], [compare_nhanes_peds()], [plot_nhanes_peds()]
#' @examples
#' sim <- sample_nhanes_peds(n = 100, ages = c(2, 8, 14), seed = 20261008)
#' head(sim)
#' attr(sim, "kernels")
#'
#' # Weight only, using every child with a measured weight
#' wt <- sample_nhanes_peds(n = 100, ages = 10, vars = "WT", seed = 1)
#' summary(wt$WT)
#' @export
sample_nhanes_peds <- function(n = 500,
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
```

Run: `Rscript -e "devtools::document()"`
Expected: `Writing 'sample_nhanes_peds.Rd'` and `NAMESPACE` gains `export(sample_nhanes_peds)`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: all tests pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 5: Commit**

```bash
git add R/nhanes_peds.R NAMESPACE man/sample_nhanes_peds.Rd man/nhanes_peds.Rd tests/testthat/test-nhanes_peds.R
git commit -m "Add sample_nhanes_peds()"
```

---

### Task 5: Simulator behaviour (fidelity, resampling, donor eligibility)

**Files:**
- Modify: `tests/testthat/test-nhanes_peds.R` (append)
- Modify: `R/nhanes_peds.R` only if a test fails

These tests pin the statistical behaviour. They should pass against the Task 4 code; if one fails, fix the implementation, not the tolerance, unless the failure shows the tolerance itself is wrong (then explain why in the commit message).

- [ ] **Step 1: Write the tests**

Append to `tests/testthat/test-nhanes_peds.R`:

```r
test_that("bandwidth_factor = 0 reproduces the donors exactly", {
  sim <- sample_nhanes_peds(n = 50, ages = 7, bandwidth_factor = 0, seed = 3)
  key_sim <- paste(sim$SOURCE_CYCLE, sim$SOURCE_SEQN)
  key_ref <- paste(nhanes_peds$CYCLE, nhanes_peds$SEQN)
  donor <- nhanes_peds[match(key_sim, key_ref), ]
  expect_equal(sim$WT, donor$WT)
  expect_equal(sim$HT, donor$HT)
})

test_that("joint simulation never uses donors with a missing measure", {
  sim <- sample_nhanes_peds(n = 500, ages = 2:17, seed = 4)
  key_sim <- paste(sim$SOURCE_CYCLE, sim$SOURCE_SEQN)
  key_ref <- paste(nhanes_peds$CYCLE, nhanes_peds$SEQN)
  donor <- nhanes_peds[match(key_sim, key_ref), ]
  expect_false(anyNA(donor$WT) || anyNA(donor$HT))
})

test_that("weight-only simulation uses more donors than joint simulation", {
  wt <- sample_nhanes_peds(n = 1, vars = "WT", seed = 5)
  both <- sample_nhanes_peds(n = 1, seed = 5)
  expect_gt(
    sum(attr(wt, "kernels")$N_NHANES),
    sum(attr(both, "kernels")$N_NHANES)
  )
})

test_that("simulated medians and correlation match the weighted reference", {
  ages <- c(2, 10, 17)
  sim <- sample_nhanes_peds(n = 4000, ages = ages, seed = 20261008)
  donors <- nhanes_donors(
    nhanes_peds, ages, c("Male", "Female"), c("WT", "HT"),
    sort(unique(nhanes_peds$CYCLE))
  )
  for (age in ages) {
    for (sx in c("Male", "Female")) {
      ref <- donors[donors$AGE == age & donors$SEX == sx, ]
      s <- sim[sim$AGE == age & sim$SEX == sx, ]
      for (v in c("WT", "HT")) {
        ref_median <- weighted_quantile(ref[[v]], ref$PROB, 0.5)
        expect_lt(abs(stats::median(s[[v]]) / ref_median - 1), 0.03)
      }
      ref_rho <- weighted_cor(log(ref$WT), log(ref$HT), ref$PROB)
      expect_lt(abs(stats::cor(log(s$WT), log(s$HT)) - ref_rho), 0.05)
    }
  }
})
```

- [ ] **Step 2: Run the tests**

Run: `Rscript -e "devtools::test(filter = 'nhanes_peds', reporter = 'summary')"`
Expected: all pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-nhanes_peds.R
git commit -m "Test sample_nhanes_peds() fidelity and donor eligibility"
```

---

### Task 6: `compare_nhanes_peds()`

**Files:**
- Create: `R/nhanes_compare.R`
- Create: `tests/testthat/test-nhanes_compare.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-nhanes_compare.R`:

```r
test_that("compare_nhanes_peds has one row per age, sex and variable", {
  sim <- sample_nhanes_peds(n = 200, ages = c(4, 11), seed = 1)
  cmp <- compare_nhanes_peds(sim)
  stats_cols <- c("N", "Mean", "SD", "P05", "Q1", "Median", "Q3", "P95")
  expect_s3_class(cmp, "tbl_df")
  expect_named(
    cmp,
    c("AGE", "SEX", "VARIABLE",
      paste0(stats_cols, "_REF"), paste0(stats_cols, "_SIM"),
      "MEAN_RATIO", "MEDIAN_RATIO", "P05_RATIO", "P95_RATIO",
      "RHO_REF", "RHO_SIM")
  )
  expect_equal(nrow(cmp), 2 * 2 * 3)
  expect_setequal(cmp$VARIABLE, c("WT", "HT", "BMI"))
  expect_true(all(cmp$N_SIM == 200))
  expect_equal(cmp$MEDIAN_RATIO, cmp$Median_SIM / cmp$Median_REF)
})

test_that("compare_nhanes_peds omits BMI and RHO for a single measure", {
  sim <- sample_nhanes_peds(n = 50, ages = 6, vars = "HT", seed = 1)
  cmp <- compare_nhanes_peds(sim)
  expect_equal(unique(cmp$VARIABLE), "HT")
  expect_false(any(c("RHO_REF", "RHO_SIM") %in% names(cmp)))
})

test_that("compare_nhanes_peds uses the simulation's cycles", {
  sim <- sample_nhanes_peds(n = 50, ages = 6, cycles = "2015-16", seed = 1)
  cmp <- compare_nhanes_peds(sim)
  n_ref <- sum(
    nhanes_peds$CYCLE == "2015-16" & nhanes_peds$AGE == 6 &
      nhanes_peds$SEX == "Male" &
      !is.na(nhanes_peds$WT) & !is.na(nhanes_peds$HT)
  )
  expect_equal(cmp$N_REF[cmp$SEX == "Male" & cmp$VARIABLE == "WT"], n_ref)
})

test_that("compare_nhanes_peds needs sample_nhanes_peds() output", {
  sim <- sample_nhanes_peds(n = 5, ages = 6, seed = 1)
  attr(sim, "cycles") <- NULL
  expect_error(compare_nhanes_peds(sim), "output of sample_nhanes_peds")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "devtools::test(filter = 'nhanes_compare', reporter = 'summary')"`
Expected: errors `could not find function "compare_nhanes_peds"`.

- [ ] **Step 3: Implement**

Create `R/nhanes_compare.R`:

```r
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
#' index by age and sex for a population from [sample_nhanes_peds()] and for
#' the NHANES reference it was drawn from. The reference uses the same
#' releases and donor eligibility as the simulation, weighted with the
#' cycle-balanced MEC weights.
#'
#' @param sim Output of [sample_nhanes_peds()], with its attributes intact.
#' @param data Reference data with the columns of [nhanes_peds]; normally the
#'   data used for the simulation.
#' @return A tibble with one row per `AGE`, `SEX` and `VARIABLE` (`WT`, `HT`
#'   and `BMI` in kg/m^2, as simulated). Columns `N`, `Mean`, `SD`, `P05`,
#'   `Q1`, `Median`, `Q3` and `P95` appear with suffix `_REF` (weighted
#'   NHANES reference) and `_SIM` (simulated), followed by the ratios
#'   `MEAN_RATIO`, `MEDIAN_RATIO`, `P05_RATIO` and `P95_RATIO` (simulated /
#'   reference). When both weight and height were simulated, `RHO_REF` and
#'   `RHO_SIM` give the correlation of log weight and log height.
#' @seealso [sample_nhanes_peds()], [plot_nhanes_peds()]
#' @examples
#' sim <- sample_nhanes_peds(n = 200, ages = c(4, 12), seed = 1)
#' compare_nhanes_peds(sim)
#' @export
compare_nhanes_peds <- function(sim, data = pmxTools::nhanes_peds) {
  cycles <- attr(sim, "cycles")
  vars <- attr(sim, "vars")
  if (is.null(cycles) || is.null(vars)) {
    stop(
      "'sim' must be the output of sample_nhanes_peds() ",
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
```

Run: `Rscript -e "devtools::document()"`
Expected: `Writing 'compare_nhanes_peds.Rd'`; `NAMESPACE` gains `export(compare_nhanes_peds)`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e "devtools::test(filter = 'nhanes_compare', reporter = 'summary')"`
Expected: all 4 tests pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 5: Commit**

```bash
git add R/nhanes_compare.R NAMESPACE man/compare_nhanes_peds.Rd tests/testthat/test-nhanes_compare.R
git commit -m "Add compare_nhanes_peds()"
```

---

### Task 7: `plot_nhanes_peds()`

**Files:**
- Modify: `R/nhanes_compare.R` (append)
- Modify: `tests/testthat/test-nhanes_compare.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-nhanes_compare.R`:

```r
test_that("plot_nhanes_peds returns a ggplot of three percentiles per source", {
  sim <- sample_nhanes_peds(n = 100, ages = c(3, 9), seed = 1)
  p <- plot_nhanes_peds(compare_nhanes_peds(sim))
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$Percentile), c("P05", "Median", "P95"))
  expect_setequal(unique(p$data$Source), c("NHANES", "Simulated"))
  expect_equal(nrow(p$data), 2 * 2 * 3 * 3 * 2)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter = 'nhanes_compare', reporter = 'summary')"`
Expected: error `could not find function "plot_nhanes_peds"`.

- [ ] **Step 3: Implement**

Append to `R/nhanes_compare.R`:

```r
#' Plot simulated against NHANES percentiles
#'
#' Plots the 5th, 50th and 95th percentiles of each simulated measure against
#' age, for the simulated population and the NHANES reference, by sex.
#'
#' @param comparison Output of [compare_nhanes_peds()].
#' @return A ggplot object.
#' @seealso [compare_nhanes_peds()], [sample_nhanes_peds()]
#' @examples
#' sim <- sample_nhanes_peds(n = 200, ages = 2:17, seed = 1)
#' plot_nhanes_peds(compare_nhanes_peds(sim))
#' @export
plot_nhanes_peds <- function(comparison) {
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
  long$VARIABLE <- factor(long$VARIABLE, levels = intersect(c("WT", "HT", "BMI"), long$VARIABLE))

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
```

(`ggplot2` is imported wholesale in `NAMESPACE`, which also provides `.data`.)

Run: `Rscript -e "devtools::document()"`
Expected: `Writing 'plot_nhanes_peds.Rd'`; `NAMESPACE` gains `export(plot_nhanes_peds)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e "devtools::test(filter = 'nhanes_compare', reporter = 'summary')"`
Expected: all 5 tests pass. Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 5: Commit**

```bash
git add R/nhanes_compare.R NAMESPACE man/plot_nhanes_peds.Rd tests/testthat/test-nhanes_compare.R
git commit -m "Add plot_nhanes_peds()"
```

---

### Task 8: NEWS, full test suite, R CMD check

**Files:**
- Modify: `NEWS.md`

- [ ] **Step 1: Add the NEWS bullet**

Under `# pmxTools (development version)` at the top of `NEWS.md`, add (with a blank line after the heading):

```markdown
* Added `sample_nhanes_peds()`, which simulates virtual pediatric populations (ages 2-17) with realistic body weight and/or height by smoothed, MEC-weighted resampling of children from four NHANES releases (2013-2023), bundled as `nhanes_peds`. `compare_nhanes_peds()` and `plot_nhanes_peds()` check the simulated population against the reference.
```

- [ ] **Step 2: Run the full test suite**

Run: `Rscript -e "devtools::test(reporter = 'summary')"`
Expected: no failures; the only skips are the two `test-plot.R` vdiffr tests ("On Windows"). Then `git checkout -- tests/testthat/_snaps`.

- [ ] **Step 3: Run R CMD check**

Run:
`Rscript -e "rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'), build_args = '--no-manual', error_on = 'never', check_dir = Sys.getenv('TEMP'))"`
Expected: 0 errors, 0 warnings. Acceptable NOTE: none expected; a NOTE about installed size or data size must be investigated (the `.rda` should be well under 1 MB). Confirm the tarball excludes `data-raw/` and `dev/`.

- [ ] **Step 4: Commit**

```bash
git add NEWS.md
git commit -m "NEWS: sample_nhanes_peds()"
```

---

## Self-review (completed while writing)

- **Spec coverage:** data (Task 1); keep records with one missing measure and
  `vars` option (Tasks 1, 3, 4, 5); 1-D/2-D bandwidth rule (Task 3);
  simulator interface, return value, attributes, seed (Task 4); validation
  errors incl. empty strata (Tasks 3, 4); fidelity and bandwidth 0 (Task 5);
  comparison incl. BMI/RHO only when both vars (Task 6); plot (Task 7); NEWS,
  `.Rbuildignore` (already committed with the spec), check (Task 8).
- **Names used consistently:** `nhanes_donors()`, `nhanes_kernel()`,
  `nhanes_sample_kernel()`, `nhanes_check_values()`, `nhanes_check_vars()`,
  `nhanes_summary()`, `weighted_quantile()`, `weighted_sd()`,
  `weighted_cor()`, `robust_scale()`; kernel list fields `y`, `p`, `H`,
  `neff`, `rho`, `cycle`, `seqn`.
