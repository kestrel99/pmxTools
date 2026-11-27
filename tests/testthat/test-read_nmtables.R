# Write a data frame as a NONMEM $TABLE file (with the "TABLE NO." banner)
local_nmtable <- function(df, dir, name) {
  path <- file.path(dir, name)
  fmt <- function(x) formatC(x, format = "E", digits = 4, width = 12)
  body <- apply(df, 1, function(row) paste(vapply(row, fmt, character(1)), collapse = ""))
  writeLines(c("TABLE NO.  1", paste(formatC(names(df), width = 12), collapse = ""), body), path)
  path
}

test_that("read_nmtables() reads a single table", {
  dir <- withr::local_tempdir()
  df <- data.frame(ID = c(1, 1, 2, 2), TIME = c(0, 1, 0, 1), DV = c(0, 5.5, 0, 4.25), PRED = c(0, 5, 0, 4))
  local_nmtable(df, dir, "sdtab1")

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  expect_equal(res, df)
})

test_that("read_nmtables() keeps the shape of a one-row table", {
  dir <- withr::local_tempdir()
  df <- data.frame(ID = 1, CL = 2.5)
  local_nmtable(df, dir, "patab1")

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  expect_equal(res, df)
})

test_that("read_nmtables() merges FIRSTONLY tables with three or more files", {
  dir <- withr::local_tempdir()
  sdtab <- data.frame(ID = c(1, 1, 2), TIME = c(0, 1, 0), DV = c(0, 5.5, 0), PRED = c(0, 5, 0))
  patab <- data.frame(ID = c(1, 1, 2), TIME = c(0, 1, 0), IPRED = c(0, 5.2, 0), CWRES = c(0, 0.3, 0))
  cotab <- data.frame(ID = c(1, 2), WT = c(70, 80), AGE = c(30, 40), SEX = c(0, 1))
  local_nmtable(sdtab, dir, "sdtab1")
  local_nmtable(patab, dir, "patab1")
  local_nmtable(cotab, dir, "cotab1")

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  expect_equal(nrow(res), 3)
  expect_equal(res$WT[order(res$ID, res$TIME)], c(70, 70, 80))
  expect_equal(res$IPRED[order(res$ID, res$TIME)], c(0, 5.2, 0))
  expect_false(anyDuplicated(names(res)) > 0)
})

test_that("read_nmtables() refuses to merge tables with no shared key", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(1, 2, 3), TIME = c(0, 1, 2)), dir, "sdtab1")
  local_nmtable(data.frame(SUBJ = c(1, 2), WT = c(70, 80)), dir, "cotab1")

  expect_error(read_nmtables(runNo = 1, directory = dir, quiet = TRUE), "ID")
})

test_that("read_nm_multi_table() keeps the shape of a one-row table", {
  path <- withr::local_tempfile()
  writeLines(c("TABLE NO.     1: First Order", " ITERATION    THETA1       THETA2", "  -1000000000  1.0E+00  2.0E+00"),
             path)

  res <- read_nm_multi_table(path)

  expect_equal(res, data.frame(ITERATION = -1e9, THETA1 = 1, THETA2 = 2))
})

test_that("read_nmtables() keeps unmatched subjects when merging and warns", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(1, 2, 3), TIME = c(0, 1, 2)), dir, "sdtab1")
  local_nmtable(data.frame(ID = c(1, 2), WT = c(70, 80)), dir, "cotab1")

  expect_warning(res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE), "ID: 3")
  expect_equal(res, data.frame(ID = c(1, 2, 3), TIME = c(0, 1, 2), WT = c(70, 80, NA)))
})

test_that("read_nmtables() merges three-column FIRSTONLY tables", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(2, 2, 1), TIME = c(0, 1, 0), DV = c(1, 2, 3), PRED = c(4, 5, 6)), dir, "sdtab1")
  local_nmtable(data.frame(ID = c(1, 2), CL = c(2.5, 3.5), V = c(20, 30)), dir, "patab1")

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  # Per-record row order is preserved
  expect_equal(res, data.frame(ID = c(2, 2, 1), TIME = c(0, 1, 0), DV = c(1, 2, 3), PRED = c(4, 5, 6),
                               CL = c(3.5, 3.5, 2.5), V = c(30, 30, 20)))
})

test_that("read_nmtables() stacks repeated tables from NSIM > 1 runs", {
  dir <- withr::local_tempdir()
  writeLines(c("TABLE NO.  1", "  ID  TIME  DV  PRED", " 1 0 1 1", " 1 1 2 2",
               "TABLE NO.  1", "  ID  TIME  DV  PRED", " 1 0 3 1", " 1 1 4 2"),
             file.path(dir, "sdtab1"))

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  expect_equal(res, data.frame(ID = 1, TIME = c(0, 1, 0, 1), DV = c(1, 2, 3, 4), PRED = c(1, 2, 1, 2)))
})

test_that("read_nmcov() labels rows and columns with the parameter names", {
  res <- read_nmcov(test_path("testData", "runR001.cov"), quiet = TRUE)
  labels <- c("THETA1", "THETA2", "THETA3", "SIGMA(1,1)", "OMEGA(1,1)", "OMEGA(2,1)", "OMEGA(2,2)")
  expect_equal(dimnames(res), list(labels, labels))
  expect_true(isSymmetric(res))
})

test_that("read_nmtables() reports missing table files only when not quiet", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(1, 2), TIME = c(0, 1)), dir, "sdtab1")

  expect_no_message(read_nmtables(runNo = 1, directory = dir, quiet = TRUE))
  msgs <- testthat::capture_messages(read_nmtables(runNo = 1, directory = dir, quiet = FALSE))
  expect_true(any(grepl("not found", msgs)))
})

test_that("read_nmtables() joins each FIRSTONLY table by ID whatever its row order", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(1, 1, 2), TIME = c(0, 1, 0), DV = c(1, 2, 3), PRED = c(1, 1, 1)), dir, "sdtab1")
  local_nmtable(data.frame(ID = c(1, 2), WT = c(70, 80), X1 = 0, X2 = 0), dir, "cotab1")
  local_nmtable(data.frame(ID = c(2, 1), AGE = c(40, 30), Y1 = 0, Y2 = 0), dir, "catab1")

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  expect_equal(res$ID, c(1, 1, 2))
  expect_equal(res$WT, c(70, 70, 80))
  expect_equal(res$AGE, c(30, 30, 40))
})

test_that("read_nmtables() refuses to bind per-record tables whose rows do not line up", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(1, 1, 2), TIME = c(0, 1, 0), DV = c(1, 2, 3), PRED = c(1, 1, 1)), dir, "sdtab1")
  local_nmtable(data.frame(ID = c(1, 2, 1), TIME = c(0, 0, 1), IPRED = c(1, 3, 2), CWRES = 0), dir, "patab1")

  expect_error(read_nmtables(runNo = 1, directory = dir, quiet = TRUE), "ID, TIME")
})

test_that("read_nmtables() compares subject IDs exactly when lining up per-record tables", {
  # 10001 and 10002 are within the relative tolerance used for rounded values
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(10001, 10002), TIME = 0, DV = c(1, 2), PRED = 1), dir, "sdtab1")
  local_nmtable(data.frame(ID = c(10002, 10001), TIME = 0, WT = c(90, 55), CWRES = 0), dir, "patab1")

  expect_error(read_nmtables(runNo = 1, directory = dir, quiet = TRUE), "shared column(s) ID,", fixed = TRUE)
})

test_that("read_nmtables() binds per-record tables written with different precision", {
  dir <- withr::local_tempdir()
  local_nmtable(data.frame(ID = c(1, 1), TIME = c(0, 1 / 3), DV = c(1, 2), PRED = c(1, 1)), dir, "sdtab1")
  writeLines(c("TABLE NO.  2", "  ID  TIME  IPRED  CWRES", " 1 0 1.5 0", " 1 0.333333333 2.5 0"),
             file.path(dir, "patab1"))

  res <- read_nmtables(runNo = 1, directory = dir, quiet = TRUE)

  expect_equal(res$IPRED, c(1.5, 2.5))
})
