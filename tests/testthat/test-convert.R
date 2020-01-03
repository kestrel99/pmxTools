library(testthat)
library(pmxTools)

context("Convert.xls PK functions")

test_that("1-compartment model", {
  t <- calc_derived_1cpt(CL=16, V=25)
  expect_equal(
    t,
    list(k10=0.64, Vss=25, thalf=1.083, alpha=0.64, trueA=0.04, fracA=1,
         V1=25, CL=16, ka=NULL, tlag=NULL)
  )
  t_vec <- calc_derived_1cpt(CL=c(16, 8), V=c(25, 25))
  expect_equal(
    t_vec,
    list(
      k10=c(0.64, 0.32),
      Vss=c(25, 25),
      thalf=c(1.083, 2.1661),
      alpha=c(0.64, 0.32),
      trueA=c(0.04, 0.04),
      fracA=1,
      V1=c(25, 25),
      CL=c(16, 8),
      ka=NULL,
      tlag=NULL
    )
  )
})

test_that("2-compartment model", {
  t <- calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
  expect_equal(
    t,
    list(
      k10=0.64, k12=0.02, k21=0.01, Vss=75, thalf_alpha=1.0497, thalf_beta=71.514, alpha=0.66031, beta=0.0096925,
      trueA=0.039981, trueB=1.8908e-05, fracA=0.99953, fracB=0.0004727,
      V1=25, V2=50, CL=16, Q2=0.5, ka=NULL, tlag=NULL
    )
  )
  t_vec <- calc_derived_2cpt(CL=c(16, 8), V1=c(25, 50), V2=c(50, 25), Q=c(0.5, 1))
  expect_equal(
    t_vec,
    list(
      k10=c(0.64, 0.16),
      k12=c(0.02, 0.02),
      k21=c(0.01, 0.04),
      Vss=c(75, 75),
      thalf_alpha=c(1.0497, 3.7367),
      thalf_beta=c(71.514, 20.090),
      alpha=c(0.66031, 0.1855),
      beta=c(0.0096925, 0.034502),
      trueA=c(0.039981, 0.019272),
      trueB=c(1.8908e-05, 7.2827e-04),
      fracA=c(0.99953, 0.96359),
      fracB=c(0.0004727, 0.0364140),
      V1=c(25, 50),
      V2=c(50, 25),
      CL=c(16, 8),
      Q2=c(0.5, 1),
      ka=NULL,
      tlag=NULL
    )
  )
})

test_that("3-compartment model", {
  t <- calc_derived_3cpt(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
  expect_equal(
    t,
    list(
      k10=1.2564, k12=11.538, k21=2.3684, k13=3.1197, k31=0.015821, Vss=4751.4, 
      thalf_alpha=0.039161, thalf_beta=1.1659, thalf_gamma=154.92,
      alpha=17.7, beta=0.59449, gamma=0.0044742,
      trueA=0.038279, trueB=0.0043467, trueC=0.0001098,
      fracA=0.89572, fracB=0.10171, fracC=0.0025692,
      V1=23.4, V2=114, V3=4614, CL=29.4, Q2=270, Q3=73, ka=NULL, tlag=NULL
    )
  )
  t_vec <-
    calc_derived_3cpt(
      CL=c(20, 10),
      V1=c(20, 10),
      V2=c(30, 40),
      V3=c(40, 100),
      Q2=c(10, 20),
      Q3=c(30, 40)
    )
  # Compares equally to:
  #t_vec1 <- calc_derived_3cpt(CL=20, V1=20, V2=30, V3=40, Q2=10, Q3=30)
  #t_vec2 <- calc_derived_3cpt(CL=10, V1=10, V2=40, V3=100, Q2=20, Q3=40)
  expect_equal(
    t_vec,
    list(
      k10=c(1, 1),
      k12=c(0.5, 2),
      k21=c(0.33333, 0.5),
      k13=c(1.5, 4),
      k31=c(0.75, 0.4),
      Vss=c(90, 150), 
      thalf_alpha=c(0.19991, 0.093988),
      thalf_beta=c(1.510, 1.484),
      thalf_gamma=c(4.4129, 11.938),
      alpha=c(3.4672, 7.3749),
      beta=c(0.45905, 0.46709),
      gamma=c(0.15707, 0.05806),
      trueA=c(0.042759, 0.094872),
      trueB=c(0.0020133, 7.8148e-05),
      trueC=c(0.0052276, 0.0050494),
      fracA=c(0.85518, 0.94872),
      fracB=c(0.040266, 0.00078148),
      fracC=c(0.10455, 0.050494),
      V1=c(20, 10),
      V2=c(30, 40),
      V3=c(40, 100),
      CL=c(20, 10),
      Q2=c(10, 20),
      Q3=c(30, 40),
      ka=NULL,
      tlag=NULL
    )
  )
})

test_that("automatic detection works", {
  expect_equal(
    expect_message(
      calc_derived(CL=16, V=25, verbose=TRUE),
      regexp="Detected 1-compartment model",
      fixed=TRUE
    ),
    calc_derived_1cpt(CL=16, V=25)
  )
  expect_equal(
    expect_message(
      calc_derived(CL=16, V1=25, V2=50, Q=0.5, verbose=TRUE),
      regexp="Detected 2-compartment model",
      fixed=TRUE
    ),
    calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
  )
  expect_equal(
    expect_message(
      calc_derived(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73, verbose=TRUE),
      regexp="Detected 3-compartment model",
      fixed=TRUE
    ),
    calc_derived_3cpt(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
  )
  expect_error(
    calc_derived(foo=1),
    regexp="Could not determine model type based on argument names.  Please check the following argument names: foo",
    fixed=TRUE
  )
})
