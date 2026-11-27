test_that("pk_curve() superposes additional doses by elapsed time on irregular grids", {
  p <- list(CL = 2.5, V = 25, ka = 2.77)
  t <- c(0, 1, 25)
  single <- function(x) calc_sd_1cmt_linear_oral_1(t = x, dose = 600, CL = 2.5, V = 25, ka = 2.77)

  res <- pk_curve(t = t, params = p, ii = 24, addl = 1)

  expect_equal(res$t, t)
  expect_equal(res$cp, c(single(0), single(1), single(25) + single(1)))
})

test_that("pk_curve() handles grids that start after zero and addl beyond the observed window", {
  p <- list(CL = 2.5, V = 25)
  t <- c(30, 49.5, 50)
  single <- function(x) calc_sd_1cmt_linear_bolus(t = x, dose = 600, CL = 2.5, V = 25)

  res <- pk_curve(t = t, model = "1cmt_bolus", params = p, ii = 24, addl = 10)

  expected <- sapply(t, function(x) {
    dose_times <- 24 * (0:10)
    sum(single(x - dose_times[dose_times <= x]))
  })
  expect_equal(res$cp, expected)
})

test_that("pk_curve() gives zero contribution from doses not yet given", {
  res <- pk_curve(t = c(0, 5, 11.99), model = "1cmt_bolus", params = list(CL = 2.5, V = 25),
                  ii = 12, addl = 3)
  single <- calc_sd_1cmt_linear_bolus(t = c(0, 5, 11.99), dose = 600, CL = 2.5, V = 25)
  expect_equal(res$cp, single)
})

test_that("pk_curve(ss = TRUE) gives the steady-state curve over the dosing interval", {
  p <- list(CL = 2.5, V1 = 25, V2 = 50, Q = 5, ka = 1)
  tad <- c(0, 0.5, 3, 11.9)

  res <- pk_curve(t = tad, model = "2cmt_oral", params = p, ii = 12, ss = TRUE)

  expect_equal(
    res$cp,
    calc_ss_2cmt_linear_oral_1(tad = tad, tau = 12, dose = 600, CL = 2.5, V1 = 25, V2 = 50, Q = 5, ka = 1)
  )
  expect_false(isTRUE(all.equal(res$cp, pk_curve(t = tad, model = "2cmt_oral", params = p, ii = 12)$cp)))
})

test_that("pk_curve(ss = TRUE) washes out after the last dose and adds addl doses", {
  p <- list(CL = 2.5, V = 25, tinf = 2)
  t <- c(0, 1, 13, 30, 40)
  ss <- function(x) calc_ss_1cmt_linear_infusion(tad = x, tau = 12, dose = 600, CL = 2.5, V = 25, tinf = 2)
  sd <- function(x) calc_sd_1cmt_linear_infusion(t = x, dose = 600, CL = 2.5, V = 25, tinf = 2)

  res <- pk_curve(t = t, model = "1cmt_infusion", params = p, ii = 12, addl = 1, ss = TRUE)

  expect_equal(res$cp, c(ss(0), ss(1), ss(13) + sd(1), ss(30) + sd(18), ss(40) + sd(28)))
})

test_that("pk_curve() at steady state matches a long run of single doses", {
  p <- list(CL = 2.5, V1 = 25, V2 = 50, V3 = 100, Q2 = 1, Q3 = 0.5)
  long <- pk_curve(t = 24 * 400 + c(0, 6, 23), model = "3cmt_bolus", params = p, ii = 24, addl = 400)
  ss <- pk_curve(t = c(0, 6, 23), model = "3cmt_bolus", params = p, ii = 24, ss = TRUE)
  expect_equal(long$cp, ss$cp, tolerance = 1e-6)
})

test_that("pk_curve() validates its inputs", {
  p <- list(CL = 2.5, V = 25, ka = 1)
  expect_error(pk_curve(t = 0:2, model = "4cmt_oral", params = p), "model")
  expect_error(pk_curve(t = 0:2, model = "1cmt_oral", params = list(CL = 2.5, V = 25)), "ka")
  expect_error(pk_curve(t = 0:2, model = "2cmt_bolus", params = list(CL = 1, V1 = 1, V2 = 1)),
               "two-compartment IV bolus")
  expect_error(pk_curve(t = c(-1, 0, 1), params = p), "t")
  expect_error(pk_curve(t = numeric(0), params = p), "t")
  expect_error(pk_curve(t = 0:2, params = p, ii = 0), "ii")
  expect_error(pk_curve(t = 0:2, params = p, addl = -1), "addl")
  expect_error(pk_curve(t = 0:2, params = p, addl = 1.5), "addl")
  expect_error(pk_curve(t = 0:2, params = p, ss = NA), "ss")
  expect_error(pk_curve(t = 0:2, params = p, dose = -1), "dose")
})
