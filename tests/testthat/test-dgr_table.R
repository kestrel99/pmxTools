test_that("dgr_table works with new API (named vector)", {
  dat <- data.frame(
    ID = 1:10,
    AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
    SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age", SEX = "Sex"))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(result$Variable[1], "N")
  expect_equal(result$Total[1], "10")
  expect_true("Age" %in% result$Variable)
  expect_true("- F" %in% result$Variable)
  expect_true("- M" %in% result$Variable)
})

test_that("dgr_table works with backward compatible API", {
  dat <- data.frame(
    ID = 1:10,
    AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
    SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- dgr_table(dat, c("AGE", "SEX"), c("Age", "Sex")),
    "deprecated"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(result$Variable[1], "N")
})

test_that("dgr_table works with grouping", {
  dat <- data.frame(
    ID = 1:10,
    AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
    SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age", SEX = "Sex"), by = "SEX")

  expect_equal(ncol(result), 4)
  expect_true("F" %in% names(result))
  expect_true("M" %in% names(result))
  expect_true("Total" %in% names(result))
  expect_equal(result$Variable[1], "N")
})

test_that("dgr_table handles continuous variables correctly", {
  dat <- data.frame(
    ID = 1:10,
    AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age"))

  age_row <- result$Total[result$Variable == "Age"]
  expect_true(grepl("[", age_row, fixed = TRUE))
  expect_true(grepl("(", age_row, fixed = TRUE))
})

test_that("dgr_table handles categorical variables correctly", {
  dat <- data.frame(
    ID = 1:10,
    SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(SEX = "Sex"))

  expect_true("- M" %in% result$Variable)
  expect_true("- F" %in% result$Variable)
  expect_true(grepl("50%", result$Total[result$Variable == "- M"]))
  expect_true(grepl("50%", result$Total[result$Variable == "- F"]))
})

test_that("dgr_table handles navars with numeric missing codes", {
  dat <- data.frame(
    ID = 1:5,
    AGE = c(-99, 52, 38, -999, 48),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age"))

  expect_true("Missing" %in% result$Variable)
  expect_true(grepl("40%", result$Total[result$Variable == "Missing"]))
})

test_that("dgr_table handles explicit numeric navars", {
  dat <- data.frame(
    ID = 1:5,
    AGE = c(-99, 52, 38, -999, 48),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age"), navars = c(-99, -999))

  expect_true("Missing" %in% result$Variable)
  expect_true(grepl("40%", result$Total[result$Variable == "Missing"]))
})

test_that("dgr_table handles grouped navars correctly", {
  dat <- data.frame(
    ID = 1:5,
    AGE = c(-99, 52, 38, -999, 48),
    SEX = c("M", "F", "M", "F", "M"),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age"), by = "SEX", navars = c(-99, -999))

  expect_true("Missing" %in% result$Variable)
})

test_that("dgr_table uses arithmetic mean when mtype = mean", {
  dat <- data.frame(
    ID = 1:20,
    AGE = 1:20,
    stringsAsFactors = FALSE
  )

  result_geom <- dgr_table(dat, c(AGE = "Age"), mtype = "geomean")
  result_arith <- dgr_table(dat, c(AGE = "Age"), mtype = "mean")

  age_row_geom <- result_geom$Total[result_geom$Variable == "Age"]
  age_row_arith <- result_arith$Total[result_arith$Variable == "Age"]

  expect_false(identical(age_row_geom, age_row_arith))
})

test_that("dgr_table respects cutoff parameter", {
  dat <- data.frame(
    ID = 1:20,
    SMALL = rep(1:2, each = 10),
    LARGE = 1:20,
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(SMALL = "Small", LARGE = "Large"), cutoff = 3)

  small_row <- result$Total[result$Variable == "Small"]
  large_row <- result$Total[result$Variable == "Large"]

  expect_false(grepl("[", small_row, fixed = TRUE))
  expect_true(grepl("[", large_row, fixed = TRUE))
})

test_that("dgr_table respects sig parameter", {
  dat <- data.frame(
    ID = 1:10,
    VAL = c(1.11111, 2.22222, 3.33333, 4.44444, 5.55555, 6.66666, 7.77777, 8.88888, 9.99999, 10.10101),
    stringsAsFactors = FALSE
  )

  result_3 <- dgr_table(dat, c(VAL = "Val"), sig = 3)
  result_5 <- dgr_table(dat, c(VAL = "Val"), sig = 5)

  val_row_3 <- result_3$Total[result_3$Variable == "Val"]
  val_row_5 <- result_5$Total[result_5$Variable == "Val"]

  expect_true(nchar(val_row_3) < nchar(val_row_5))
})

test_that("dgr_table handles single row data", {
  dat <- data.frame(
    ID = 1,
    AGE = 45,
    SEX = "M",
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age", SEX = "Sex"))

  expect_equal(result$Total[1], "1")
  expect_true("Age" %in% result$Variable)
  expect_true("- M" %in% result$Variable)
})

test_that("dgr_table handles all same values", {
  dat <- data.frame(
    ID = 1:5,
    VAL = rep(100, 5),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(VAL = "Value"))

  expect_true("Value" %in% result$Variable)
  expect_true("- 100" %in% result$Variable)
  expect_true(grepl("100", result$Total[result$Variable == "- 100"]))
})

test_that("dgr_table handles grouped continuous variable", {
  dat <- data.frame(
    ID = 1:10,
    AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
    SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(AGE = "Age"), by = "SEX")

  expect_equal(ncol(result), 4)
  age_row <- result$Total[result$Variable == "Age"]
  expect_true(grepl("[", age_row, fixed = TRUE))
  expect_true(grepl("(", age_row, fixed = TRUE))
})

test_that("dgr_table handles grouped categorical variable", {
  dat <- data.frame(
    ID = 1:10,
    SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F"),
    GRP = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(SEX = "Sex"), by = "GRP")

  expect_equal(ncol(result), 4)
  expect_true("A" %in% names(result))
  expect_true("B" %in% names(result))
})

test_that("dgr_table handles empty groups", {
  dat <- data.frame(
    ID = 1:5,
    VAL = c(10, 20, 30, 40, 50),
    GRP = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )

  result <- dgr_table(dat, c(VAL = "Value"), by = "GRP")

  expect_equal(ncol(result), 4)
  expect_true("A" %in% names(result))
  expect_true("B" %in% names(result))
})
