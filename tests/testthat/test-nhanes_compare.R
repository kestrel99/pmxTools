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

test_that("plot_nhanes_peds returns a ggplot of three percentiles per source", {
  sim <- sample_nhanes_peds(n = 100, ages = c(3, 9), seed = 1)
  p <- plot_nhanes_peds(compare_nhanes_peds(sim))
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$Percentile), c("P05", "Median", "P95"))
  expect_setequal(unique(p$data$Source), c("NHANES", "Simulated"))
  expect_equal(nrow(p$data), 2 * 2 * 3 * 3 * 2)
})
