library(testthat)
library(pmxTools)

context("Closed-form solutions")

test_that("1-compartment linear, IV bolus, single-dose", {
  t <- calc_sd_1cmt_linear_bolus(t=seq(0, 24, by=3), CL=6, V=25, dose=600)
  expect_equal(signif(t, 4), c(24.0, 11.68, 5.686, 2.768, 1.347, 0.6558, 0.3192, 0.1554, 0.07563))
})

test_that("1-compartment linear, infusion, single-dose", {
  t <- calc_sd_1cmt_linear_infusion(t=seq(0, 24, by=3), CL=6, V=25, dose=600, tinf=1)
  expect_equal(signif(t, 4), c(0, 13.2, 6.427, 3.128, 1.523, 0.7412, 0.3608, 0.1756, 0.08547))
})

test_that("1-compartment linear, zero-order oral, single-dose", {
  t <- calc_sd_1cmt_linear_oral_0(t=seq(0, 24, by=3), CL=6, V=25, dur=1.5, dose=600)
  expect_equal(signif(t, 4), c(0, 14.06, 6.845, 3.332, 1.622, 0.7893, 0.3842, 0.187, 0.09103))
})

test_that("1-compartment linear, first-order oral, single-dose", {
  t <- calc_sd_1cmt_linear_oral_1(t=seq(0, 24, by=3), CL=6, V=25, ka=1.1, dose=600)
  expect_equal(signif(t, 4), c(0, 13.81, 7.231, 3.53900, 1.723, 0.8388, 0.4083, 0.1987, 0.09673))
})

test_that("1-compartment linear, first-order oral with lag time, single-dose", {
  t <- calc_sd_1cmt_linear_oral_1_lag(t=seq(0, 24, by=3), CL=6, V=25, ka=1.1, dose=600, tlag=2)
  expect_equal(signif(t, 4), c(0, 13.93, 11.38, 5.707, 2.784, 1.356, 0.6598, 0.3212, 0.1563))
})

test_that("2-compartment linear, infusion, single-dose", {
  t <- calc_sd_2cmt_linear_infusion(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,dose = 10, tinf = 1)
  expect_equal(signif(t, 4), c(0, 0.1855, 0.05658, 0.01783, 0.006158, 0.002617, 0.001519, 0.001156, 0.001015))
})

test_that("1-compartment linear, IV bolus, steady-state", {
  t <- calc_ss_1cmt_linear_bolus(tad=seq(0, 24, by=3), CL=6, V=25, dose=600, tau=24)
  expect_equal(signif(t, 4), c(24.08, 11.72, 5.704, 2.777, 1.351, 0.6578, 0.3202, 0.1559, 0.07587))
})

test_that("1-compartment linear, infusion, steady-state", {
  t <- calc_ss_1cmt_linear_infusion(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, tinf=1, tau=24)
  expect_equal(signif(t, 4), c(4.292, 23.03, 18.12, 14.25, 11.21, 8.819, 6.937, 5.457, 4.292))
})

test_that("1-compartment linear, zero-order oral, steady-state", {
  t <- calc_ss_1cmt_linear_oral_0(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, dur=1, tau=24)
  expect_equal(signif(t, 4), c(4.292, 23.03, 18.12, 14.25, 11.21, 8.819, 6.937, 5.457, 4.292))
})

test_that("1-compartment linear, first-order oral, steady-state", {
  t <- calc_ss_1cmt_linear_oral_1(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, ka=0.25, tau=24)
  expect_equal(signif(t, 4), c(5.976, 15.82, 17.7, 16.4, 14.07, 11.62, 9.406, 7.522, 5.976))
})

test_that("1-compartment linear, first-order oral with lag time, steady-state", {
  t <- calc_ss_1cmt_linear_oral_1_lag(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, ka=0.25, tlag=0.75, tau=24)
  expect_equal(signif(t, 4), c(6.332, 14.38, 17.65, 16.88, 14.69, 12.22, 9.931, 7.961, 6.332))
})

test_that("2-compartment linear, infusion, steady-state", {
  t <- calc_ss_2cmt_linear_infusion(tad = seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, tinf = 1, tau = 12)
  expect_equal(signif(t, 4), c(0.01191, 0.1936, 0.0633, 0.02395, 0.01191, 0.008084, 0.006728, 0.006125, 0.005756))
})
