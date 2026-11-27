# A forward-only PsN SCM log with one step and two candidate relationships
local_scm_dir <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  writeLines(c(
    "Model directory C:/scm/m1",
    "MODEL            TEST     BASE OFV     NEW OFV    TEST OFV (DROP)    GOAL     dDF    SIGNIFICANT PVAL",
    "CLWT-2           PVAL    1000.000   990.000            10.000  >   3.841    1  YES!  0.001565",
    "CLAGE-2          PVAL    1000.000   999.000             1.000  >   3.841    1        0.317311"
  ), file.path(dir, "scmlog.txt"))
  writeLines("", file.path(dir, "short_scmlog.txt"))
  dir
}

test_that("read_scm() parses a forward-only log", {
  res <- read_scm(local_scm_dir())

  expect_named(res, c("forward", "forwardSummary", "forwardP"))
  expect_equal(res$forward$Model, c("CLWT-2", "CLAGE-2"))
  expect_equal(res$forward$Drop, c(10, 1))
  expect_equal(res$forward$Significant, c("1", "0"))
  expect_equal(res$forwardSummary$BestModel, "CLWT-2")
  expect_equal(res$forwardP, signif(1 - stats::pchisq(3.841, df = 1), 2))
})

test_that("read_scm() ignores same-named objects in the global environment", {
  assign("scmb", data.frame(junk = 1), envir = globalenv())
  withr::defer(rm("scmb", envir = globalenv()))

  res <- read_scm(local_scm_dir())

  expect_null(res$backward)
  expect_named(res, c("forward", "forwardSummary", "forwardP"))
})

test_that("plot_scm() output does not depend on same-named global objects", {
  skip_if_not_installed("DiagrammeR")
  dir <- withr::local_tempdir()
  writeLines(c(
    "Model directory C:/scm/m1",
    "MODEL            TEST     BASE OFV     NEW OFV    TEST OFV (DROP)    GOAL     dDF    SIGNIFICANT PVAL",
    "CLWT-2           PVAL    1000.000   990.000            10.000  >   3.841    1  YES!  0.001565",
    "CLAGE-2          PVAL    1000.000   999.000             1.000  >   3.841    1        0.317311",
    "Model directory C:/scm/m2",
    "MODEL            TEST     BASE OFV     NEW OFV    TEST OFV (DROP)    GOAL     dDF    SIGNIFICANT PVAL",
    "CLAGE-2          PVAL     990.000   985.000             5.000  >   3.841    1  YES!  0.025347"
  ), file.path(dir, "scmlog.txt"))
  writeLines("", file.path(dir, "short_scmlog.txt"))

  clean <- plot_scm(dir)
  assign("scm_bwd", "junk", envir = globalenv())
  withr::defer(rm("scm_bwd", envir = globalenv()))
  polluted <- plot_scm(dir)

  expect_s3_class(clean, "grViz")
  expect_equal(polluted$x$diagram, clean$x$diagram)
})

test_that("read_scm() records each step's model directory, whatever the path style", {
  dir <- withr::local_tempdir()
  writeLines(c(
    "Model directory /home/user/scm_run/m1",
    "MODEL            TEST     BASE OFV     NEW OFV    TEST OFV (DROP)    GOAL     dDF    SIGNIFICANT PVAL",
    "CLWT-2           PVAL    1000.000   990.000            10.000  >   3.841    1  YES!  0.001565",
    "Model directory d:\\scm_run\\m2",
    "MODEL            TEST     BASE OFV     NEW OFV    TEST OFV (DROP)    GOAL     dDF    SIGNIFICANT PVAL",
    "CLAGE-2          PVAL     990.000   985.000             5.000  >   3.841    1  YES!  0.025347"
  ), file.path(dir, "scmlog.txt"))
  writeLines("", file.path(dir, "short_scmlog.txt"))

  res <- read_scm(dir)

  expect_equal(res$forward$Dir, c("/home/user/scm_run/m1", "d:\\scm_run\\m2"))
  expect_equal(as.character(res$forward$Step), c("1", "2"))
})

test_that("plot_scm() draws an SCM that stopped after one forward step", {
  skip_if_not_installed("DiagrammeR")

  res <- plot_scm(local_scm_dir())

  expect_s3_class(res, "grViz")
  expect_match(res$x$diagram, "CLWT-2")
  expect_match(res$x$diagram, "CLAGE-2")
})
