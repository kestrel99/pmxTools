test_that("cut_quantile works with single variable", {
  dat <- data.frame(
    ID = 1:100,
    AGE = rnorm(100, 50, 10),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, "AGE", n_groups = 4)

  expect_true("AGEQ4Q" %in% names(result))
  expect_true("AGEQ4C" %in% names(result))
  expect_equal(nrow(result), 100)
  expect_false(any(is.na(result$AGEQ4Q)))
})

test_that("cut_quantile works with multiple cuts", {
  dat <- data.frame(
    ID = 1:100,
    AGE = rnorm(100, 50, 10),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, "AGE", n_groups = c(4, 3))

  expect_true("AGEQ4Q" %in% names(result))
  expect_true("AGET3Q" %in% names(result))
  expect_true("AGEQ4C" %in% names(result))
  expect_true("AGET3C" %in% names(result))
})

test_that("cut_quantile works with multiple variables", {
  dat <- data.frame(
    ID = 1:100,
    AGE = rnorm(100, 50, 10),
    WT = rnorm(100, 70, 10),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, c("AGE", "WT"), n_groups = 4)

  expect_true("AGEQ4Q" %in% names(result))
  expect_true("WTQ4Q" %in% names(result))
})

test_that("cut_quantile handles missing codes", {
  dat <- data.frame(
    ID = 1:10,
    VAL = c(1:5, -99, -999, 8:10),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, "VAL", n_groups = 4, missing_codes = c(-99, -999))

  expect_equal(sum(is.na(result$VALQ4Q)), 2)
  expect_equal(sum(result$VALQ4C == "BLQ"), 2)
})

test_that("cut_quantile skips zero-range bins with warning", {
  dat <- data.frame(
    ID = 1:10,
    VAL = c(rep(5, 5), 10, 15, 20, 25, 30),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, "VAL", n_groups = 4)
  # Should skip the cut due to zero-range bins
  # The original column should still be there but no binning columns added
  expect_false("VALQ4Q" %in% names(result))
})

test_that("cut_quantile handles verbose output", {
  dat <- data.frame(
    ID = 1:100,
    AGE = rnorm(100, 50, 10),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, "AGE", n_groups = 4, verbose = TRUE)

  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("summary" %in% names(result))
})

test_that("cut_quantile with id uses subject-level data", {
  dat <- data.frame(
    ID = c(1, 1, 2, 2, 3, 3),
    SUBJID = c(1, 1, 2, 2, 3, 3),
    AGE = c(50, 50, 60, 60, 70, 70),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, "AGE", n_groups = 4, id = "SUBJID")

  expect_equal(nrow(result), 6)
  expect_true("AGEQ4Q" %in% names(result))
})

test_that("cut_quantile stops with conflicting values", {
  dat <- data.frame(
    ID = c(1, 1, 2, 2),
    SUBJID = c(1, 1, 2, 2),
    AGE = c(50, 55, 60, 60),
    stringsAsFactors = FALSE
  )

  expect_error(
    cut_quantile(dat, "AGE", n_groups = 4, id = "SUBJID"),
    "conflicting values"
  )
})

test_that("cut_quantile handles named list input", {
  dat <- data.frame(
    ID = 1:100,
    AGE = rnorm(100, 50, 10),
    WT = rnorm(100, 70, 10),
    stringsAsFactors = FALSE
  )

  result <- cut_quantile(dat, list(AGE = c(4, 3), WT = 4))

  expect_true("AGEQ4Q" %in% names(result))
  expect_true("AGET3Q" %in% names(result))
  expect_true("WTQ4Q" %in% names(result))
})
