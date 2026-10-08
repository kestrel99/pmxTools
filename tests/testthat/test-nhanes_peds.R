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
