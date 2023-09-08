library(testthat)
library(pmxTools)

### 1-compartment

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

test_that("1-compartment linear, zero-order oral with lag time, single-dose", {
  t <- calc_sd_1cmt_linear_oral_0_lag(t=seq(0, 24, by=3), CL=6, V=25, dur=1.5, dose=600, tlag = 0.8)
  expect_equal(signif(t, 4), c(0, 17.0400, 8.2930, 4.0370, 1.9650, 0.9564, 0.4655, 0.2266, 0.1103))
})

test_that("1-compartment linear, first-order oral, single-dose", {
  t <- calc_sd_1cmt_linear_oral_1(t=seq(0, 24, by=3), CL=6, V=25, ka=1.1, dose=600)
  expect_equal(signif(t, 4), c(0, 13.81, 7.231, 3.53900, 1.723, 0.8388, 0.4083, 0.1987, 0.09673))
})

test_that("1-compartment linear, first-order oral with lag time, single-dose", {
  t <- calc_sd_1cmt_linear_oral_1_lag(t=seq(0, 24, by=3), CL=6, V=25, ka=1.1, dose=600, tlag=2)
  expect_equal(signif(t, 4), c(0, 13.93, 11.38, 5.707, 2.784, 1.356, 0.6598, 0.3212, 0.1563))
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

test_that("1-compartment linear, zero-order oral with lag time, steady-state", {
  t <- calc_ss_1cmt_linear_oral_0_lag(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, dur=1, tau=24, tlag=0.8)
  expect_equal(signif(t, 4), c(4.576, 24.550, 19.310, 15.190, 11.950, 9.401, 7.395, 5.817, 4.576))
})

test_that("1-compartment linear, first-order oral, steady-state", {
  t <- calc_ss_1cmt_linear_oral_1(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, ka=0.25, tau=24)
  expect_equal(signif(t, 4), c(5.976, 15.82, 17.7, 16.4, 14.07, 11.62, 9.406, 7.522, 5.976))
})

test_that("1-compartment linear, first-order oral with lag time, steady-state", {
  t <- calc_ss_1cmt_linear_oral_1_lag(tad=seq(0, 24, by=3), CL=2, V=25, dose=600, ka=0.25, tlag=0.75, tau=24)
  expect_equal(signif(t, 4), c(6.332, 14.38, 17.65, 16.88, 14.69, 12.22, 9.931, 7.961, 6.332))
})

### 2-compartment

test_that("2-compartment linear, IV bolus, single-dose", {
  t <- calc_sd_2cmt_linear_bolus(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10)
  expect_equal(signif(t, 4), c(0.5, 0.151000, 0.046220, 0.014710, 0.005216, 0.002329, 0.001427, 0.001123, 0.001000))
})

test_that("2-compartment linear, infusion, single-dose", {
  t <- calc_sd_2cmt_linear_infusion(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, tinf = 1)
  expect_equal(signif(t, 4), c(0, 0.1855, 0.05658, 0.01783, 0.006158, 0.002617, 0.001519, 0.001156, 0.001015))
})

test_that("2-compartment linear, zero-order oral, single-dose", {
  t <- calc_sd_2cmt_linear_oral_0(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, dur = 0.75)
  expect_equal(signif(t, 4), c(0.000000, 0.176000, 0.053720, 0.016970, 0.005898, 0.002538, 0.001494, 0.001147, 0.001011))
})

test_that("2-compartment linear, zero-order oral with lag time, single-dose", {
  t <- calc_sd_2cmt_linear_oral_0_lag(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, dur = 0.75, tlag = 1.2)
  expect_equal(signif(t, 4), c(0.000000, 0.284000, 0.086150, 0.026720, 0.008842, 0.003436, 0.001778, 0.001246, 0.001053))
})

test_that("2-compartment linear, first-order oral, single-dose", {
  t <- calc_sd_2cmt_linear_oral_1(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, ka = 1.5)
  expect_equal(signif(t, 4), c(0.000000, 0.198100, 0.062550, 0.019650, 0.006707, 0.002784, 0.001571, 0.001174, 0.001022))
})

test_that("2-compartment linear, first-order oral with lag time, single-dose", {
  t <- calc_sd_2cmt_linear_oral_1_lag(t=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, ka = 1.5, tlag = 1.2)
  expect_equal(signif(t, 4), c(0.000000, 0.286200, 0.100100, 0.031050, 0.010150, 0.003831, 0.001900, 0.001286, 0.001068))
})


test_that("2-compartment linear, IV bolus, steady-state", {
  t <- calc_ss_2cmt_linear_bolus(tad=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, tau = 12)
  expect_equal(signif(t, 4), c(0.510900, 0.158800, 0.052810, 0.020770, 0.010920, 0.007751, 0.006595, 0.006053, 0.005704))
})

test_that("2-compartment linear, infusion, steady-state", {
  t <- calc_ss_2cmt_linear_infusion(tad = seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, tinf = 1, tau = 12)
  expect_equal(signif(t, 4), c(0.01191, 0.1936, 0.0633, 0.02395, 0.01191, 0.008084, 0.006728, 0.006125, 0.005756))
})

test_that("2-compartment linear, zero-order oral, steady-state", {
  t <- calc_ss_2cmt_linear_oral_0(tad=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, dur = 0.75, tau=6)
  expect_equal(signif(t, 4), c(0.07205, 0.20710, 0.07205, 0.03107, 0.01833, 0.01410, 0.01244, 0.01156, 0.01094))
})

test_that("2-compartment linear, zero-order oral with lag time, steady-state", {
  t <- calc_ss_2cmt_linear_oral_0_lag(tad=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, dur = 0.75, tlag = 1.2, tau = 6)
  expect_equal(signif(t, 4), c(0.10790, 0.32600, 0.10790, 0.04202, 0.02179, 0.01530, 0.01295, 0.01186, 0.01117))
})

test_that("2-compartment linear, first-order oral, steady-state", {
  t <- calc_ss_2cmt_linear_oral_1(tad=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, ka = 1.5, tau=6)
  expect_equal(signif(t, 4), c(0.08183, 0.23220, 0.08183, 0.03407, 0.01928, 0.01442, 0.01257, 0.01164, 0.01100))
})

test_that("2-compartment linear, first-order oral with lag time, steady-state", {
  t <- calc_ss_2cmt_linear_oral_1_lag(tad=seq(0, 24, by=3), CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10, ka = 1.5, tlag = 1.2, tau=6)
  expect_equal(signif(t, 4), c(0.12340, 0.33300, 0.12340, 0.04684, 0.02328, 0.01579, 0.01314, 0.01196, 0.01124))
})


### 3-compartment

test_that("3-compartment linear, IV bolus, single-dose", {
  t <- calc_sd_3cmt_linear_bolus(t=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280)
  expect_equal(signif(t, 4), c(13.93, 0.14400, 0.07859, 0.05187, 0.03945, 0.03250, 0.02782, 0.02421, 0.02122))
})

test_that("3-compartment linear, infusion, single-dose", {
  t <- calc_sd_3cmt_linear_infusion(t=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, tinf = 1)
  expect_equal(signif(t, 4), c(0, 0.16290, 0.08594, 0.05502, 0.04103, 0.03347, 0.02851, 0.02476, 0.02168))
})

test_that("3-compartment linear, zero-order oral, single-dose", {
  t <- calc_sd_3cmt_linear_oral_0(t=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, dur = 0.75)
  expect_equal(signif(t, 4), c(0.000000, 0.15780, 0.08396, 0.05418, 0.04062, 0.03322, 0.02834, 0.02462, 0.02157))
})

test_that("3-compartment linear, zero-order oral with lag time, single-dose", {
  t <- calc_sd_3cmt_linear_oral_0_lag(t=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, dur = 0.75, tlag = 1.2)
  expect_equal(signif(t, 4), c(0.000000, 0.21390, 0.10550, 0.06318, 0.04495, 0.03575, 0.03010, 0.02601, 0.02273))
})

test_that("3-compartment linear, first-order oral, single-dose", {
  expect_warning(t <- calc_sd_3cmt_linear_oral_1(t=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, ka = 1.5))
  expect_equal(signif(t, 4), c(0.000000, 0.18840, 0.08987, 0.05655, 0.04174, 0.03386, 0.02878, 0.02497, 0.02185))
})

test_that("3-compartment linear, first-order oral with lag time, single-dose", {
  expect_warning(t <- calc_sd_3cmt_linear_oral_1_lag(t=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, ka = 1.5, tlag = 1.2))
  expect_equal(signif(t, 4), c(0.000000, 0.33090, 0.11500, 0.06653, 0.04644, 0.03653, 0.03060, 0.02639, 0.02303))
})

test_that("3-compartment linear, IV bolus, steady-state", {
  t <- calc_ss_3cmt_linear_bolus(tad=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, tau=12)
  expect_equal(signif(t, 4), c(14.02000, 0.22330, 0.14750, 0.11230, 0.09257, 0.07924, 0.06896, 0.06044, 0.05312))
})

test_that("3-compartment linear, infusion, steady-state", {
  t <- calc_ss_3cmt_linear_infusion(tad=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, tinf = 1, tau=12)
  expect_equal(signif(t, 4), c(0.09530, 0.24410, 0.15650, 0.11680, 0.09530, 0.08121, 0.07054, 0.06177, 0.05427))
})

test_that("3-compartment linear, zero-order oral, steady-state", {
  t <- calc_ss_3cmt_linear_oral_0(tad=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, dur = 0.75, tau=12)
  expect_equal(signif(t, 4), c(0.09459, 0.23850, 0.15410, 0.11560, 0.09459, 0.08071, 0.07014, 0.06143, 0.05398))
})

test_that("3-compartment linear, zero-order oral with lag time, steady-state", {
  t <- calc_ss_3cmt_linear_oral_0_lag(tad=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, dur = 0.75, tlag = 1.2, tau=12)
  expect_equal(signif(t, 4), c(0.10180, 0.29960, 0.17960, 0.12790, 0.10180, 0.08573, 0.07409, 0.06474, 0.05683))
})

test_that("3-compartment linear, first-order oral, steady-state", {
  expect_warning(t <- calc_ss_3cmt_linear_oral_1(tad=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, ka = 1.5, tau=12))
  expect_equal(signif(t, 4), c(0.09642, 0.27040, 0.16100, 0.11880, 0.09642, 0.08196, 0.07112, 0.06225, 0.05468))
})

test_that("3-compartment linear, first-order oral with lag time, steady-state", {
  expect_warning(t <- calc_ss_3cmt_linear_oral_1_lag(tad=seq(0, 24, by=3), CL = 87.6, V1 = 20.1, V2 = 186, V3=749, Q2 = 111, Q3 = 53.4, dose = 280, ka = 1.5, tlag = 1.2, tau=12))
  expect_equal(signif(t, 4), c(0.10400, 0.41810, 0.19010, 0.13210, 0.10400, 0.08716, 0.07516, 0.06561, 0.05757))
})

### PK curves


test_that("PK curves", {

  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="1cmt_bolus", ii=12, addl=5,
                        params=list(CL=2.5, V=25)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="1cmt_oral", ii=12, addl=5,
                        params=list(CL=2.5, V=25, ka=1)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="1cmt_infusion", ii=12, addl=5,
                        params=list(CL=2.5, V=25, tinf=3)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="2cmt_bolus", ii=24, addl=5,
                        params=list(CL=2.5, V1=25, V2=50, Q=5)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="2cmt_oral", ii=24, addl=5,
                        params=list(CL=2.5, V1=25, V2=50, Q=5, ka=1)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="2cmt_infusion", ii=24, addl=5,
                        params=list(CL=2.5, V1=25, V2=50, Q=5, tinf=3)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="3cmt_bolus", ii=24, addl=5,
                        params=list(CL=0.25, V1=25, V2=50, V3=100, Q2=1, Q3=0.1)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="3cmt_oral", ii=24, addl=5,
                        params=list(CL=0.25, V1=25, V2=50, V3=100, Q2=1, Q3=0.1, ka=1)))
  expect_snapshot(pk_curve(t=seq(0,72,by=0.1), model="3cmt_infusion", ii=24, addl=5,
                        params=list(CL=0.25, V1=25, V2=50, V3=100, Q2=1, Q3=0.1, tinf=3)))
})

test_that("calc_derived_xxx works with vector inputs (#29)", {
  mydata <-
    data.frame(
      CL = c(0.25, 0.5),
      Vcentral = c(1, 1),
      Vperiph1 = c(0.25, 0.25),
      Vperiph2 = c(0.25, 0.25),
      Qcp1 = c(0.25, 0.25),
      Qcp2 = c(0.25, 0.25),
      Ka = c(1, 1)
    )

  expect_type(
    calc_derived_1cpt(
      CL=mydata$CL, V1=mydata$Vcentral,
      ka=mydata$Ka,
      sigdig=Inf
    ),
    "list"
  )
  expect_type(
    calc_derived_2cpt(
      CL=mydata$CL, V1=mydata$Vcentral,
      V2=mydata$Vperiph1, Q2=mydata$Qcp1,
      ka=mydata$Ka,
      sigdig=Inf
    ),
    "list"
  )
  expect_type(
    calc_derived_3cpt(
      CL=mydata$CL, V1=mydata$Vcentral,
      V2=mydata$Vperiph1, Q2=mydata$Qcp1,
      V3=mydata$Vperiph2, Q3=mydata$Qcp2,
      ka=mydata$Ka,
      sigdig=Inf
    ),
    "list"
  )
})
