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

test_that("nhanes_donors gives each cycle an equal share within a stratum", {
  donors <- nhanes_donors(
    nhanes_peds, ages = 8, sex = "Female", vars = "WT",
    cycles = unique(nhanes_peds$CYCLE)
  )
  shares <- tapply(donors$PROB, donors$CYCLE, sum)
  expect_equal(as.vector(shares), rep(0.25, 4))
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
  expect_error(sample_nhanes_peds(ages = c(1, 18)), "'ages' value(s) not in the reference data: 1, 18", fixed = TRUE)
  expect_error(sample_nhanes_peds(sex = "M"), "'sex' value(s) not in the reference data: M", fixed = TRUE)
  expect_error(sample_nhanes_peds(cycles = "1999-00"), "'cycles' value(s) not in the reference data: 1999-00", fixed = TRUE)
  expect_error(
    sample_nhanes_peds(data = nhanes_peds[c("AGE", "SEX")]),
    "'data' is missing column(s): CYCLE, SEQN, WT, HT, MEC_WT",
    fixed = TRUE
  )
})

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
