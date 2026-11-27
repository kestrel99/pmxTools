# Copy the NONMEM XML fixture with its lines transformed by `edit`
local_nm_xml <- function(edit = identity, env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".xml", .local_envir = env)
  writeLines(edit(readLines(test_path("testData", "runR001.xml"))), path)
  path
}

test_that("read_nm() accepts valid XML regardless of trailing whitespace or line endings", {
  expected <- read_nm(test_path("testData", "runR001.xml"), quiet = TRUE)

  expect_equal(read_nm(local_nm_xml(function(x) c(x, "", "   ")), quiet = TRUE), expected)
  lf <- local_nm_xml()
  writeBin(charToRaw(paste0(paste(readLines(lf), collapse = "\n"), "\n")), lf)
  expect_equal(read_nm(lf, quiet = TRUE), expected)
})

test_that("read_nm() warns and returns NULL for files that are not valid NONMEM XML", {
  expect_warning(res <- read_nm(local_nm_xml(function(x) c("  ", x)), quiet = TRUE), "Invalid XML")
  expect_null(res)

  truncated <- local_nm_xml(function(x) x[seq_len(length(x) - 5)])
  expect_warning(res <- read_nm(truncated, quiet = TRUE), "Invalid XML")
  expect_null(res)

  empty <- withr::local_tempfile(fileext = ".xml")
  file.create(empty)
  expect_warning(res <- read_nm(empty, quiet = TRUE), "Invalid XML")
  expect_null(res)

  other_root <- withr::local_tempfile(fileext = ".xml")
  writeLines(c('<?xml version="1.0"?>', "<results><a>1</a></results>"), other_root)
  expect_warning(res <- read_nm(other_root, quiet = TRUE), "nm:output")
  expect_null(res)
})
