# Deterministic test data. Reference values come from polycor 0.8-2.

ordinal_pair <- function() {
  tab <- matrix(c(30, 15, 5, 10, 25, 15, 5, 10, 35), 3)
  list(
    x = factor(rep(rep(c("a", "b", "c"), 3), times = tab),
               levels = c("a", "b", "c"), ordered = TRUE),
    y = factor(rep(rep(c("lo", "mid", "hi"), each = 3), times = tab),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
}

serial_pair <- function() {
  n <- 200
  u <- stats::qnorm(stats::ppoints(n))
  list(
    x = u,
    y = cut(u + 0.8 * sin(1:n), c(-Inf, -0.5, 0.5, Inf),
            labels = c("L", "M", "H"), ordered_result = TRUE)
  )
}

mixed_data <- function() {
  o <- ordinal_pair()
  data.frame(
    num1 = seq_along(o$x) + sin(seq_along(o$x)),
    num2 = cos(seq_along(o$x)),
    ord1 = o$x,
    ord2 = o$y
  )
}

test_that("polychoric matches the polycor two-step estimate", {
  o <- ordinal_pair()
  expect_equal(polychoric(o$x, o$y), 0.667189, tolerance = 1e-5)
})

test_that("polyserial matches the polycor ML estimate", {
  s <- serial_pair()
  expect_equal(polyserial(s$x, s$y), 0.870706, tolerance = 1e-3)
})

test_that("polychoric and polyserial recover a known latent correlation", {
  # n is large so sampling error (SD about 0.01) stays well inside 0.03
  set.seed(1)
  n <- 20000
  z1 <- stats::rnorm(n)
  z2 <- 0.5 * z1 + sqrt(1 - 0.5^2) * stats::rnorm(n)
  o1 <- cut(z1, c(-Inf, -0.5, 0.5, Inf), ordered_result = TRUE)
  o2 <- cut(z2, c(-Inf, 0, 1, Inf), ordered_result = TRUE)
  expect_lt(abs(polychoric(o1, o2) - 0.5), 0.03)
  expect_lt(abs(polyserial(z1, o2) - 0.5), 0.03)
})

test_that("polychoric and polyserial drop missing values", {
  o <- ordinal_pair()
  x <- o$x
  x[1:3] <- NA
  expect_equal(polychoric(x, o$y), polychoric(o$x[-(1:3)], o$y[-(1:3)]))

  s <- serial_pair()
  num <- s$x
  num[5] <- NA
  expect_equal(polyserial(num, s$y), polyserial(s$x[-5], s$y[-5]))
})

test_that("degenerate pairs give NA", {
  o <- ordinal_pair()
  one_level <- factor(rep("a", length(o$y)), ordered = TRUE)
  expect_identical(polychoric(one_level, o$y), NA_real_)
  expect_identical(polyserial(rep(1, length(o$y)), o$y), NA_real_)
})

test_that("pair functions check their inputs", {
  o <- ordinal_pair()
  expect_error(polychoric(as.integer(o$x), o$y), "'x' and 'y' must be factors")
  expect_error(polyserial(o$x, o$y), "'x' must be numeric")
  expect_error(polyserial(as.numeric(o$x), 1:150), "'y' must be a factor")
  expect_error(polychoric(o$x, o$y[-1]), "same length")
})

test_that("pair functions warn about unordered factors", {
  o <- ordinal_pair()
  expect_warning(
    polychoric(factor(o$x, ordered = FALSE), o$y),
    "unordered factor 'x' as ordinal",
    fixed = TRUE
  )
  s <- serial_pair()
  expect_warning(
    polyserial(s$x, factor(s$y, ordered = FALSE)),
    "unordered factor 'y' as ordinal",
    fixed = TRUE
  )
})

test_that("cor_types classifies every pair", {
  types <- cor_types(mixed_data())
  expect_equal(types["num1", "num2"], "Pearson")
  expect_equal(types["num1", "ord1"], "Polyserial")
  expect_equal(types["ord1", "num2"], "Polyserial")
  expect_equal(types["ord1", "ord2"], "Polychoric")
  expect_equal(unname(diag(types)), c("Pearson", "Pearson", "Polychoric", "Polychoric"))
})

test_that("het_cor combines the three estimators", {
  d <- mixed_data()
  R <- het_cor(d)
  expect_equal(dim(R), c(4L, 4L))
  expect_equal(dimnames(R), list(names(d), names(d)))
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_equal(R, t(R))
  expect_equal(R["num1", "num2"], stats::cor(d$num1, d$num2))
  expect_equal(R["ord1", "ord2"], polychoric(d$ord1, d$ord2))
  expect_equal(R["num1", "ord1"], polyserial(d$num1, d$ord1))
})

test_that("het_cor equals cor() for all-numeric data", {
  d <- mixed_data()[c("num1", "num2")]
  d$num3 <- d$num1 * d$num2
  expect_equal(het_cor(d), stats::cor(d))
})

test_that("het_cor uses pairwise-complete cases", {
  d <- mixed_data()
  d$num1[1:10] <- NA
  d$ord2[11:20] <- NA
  R <- het_cor(d)
  ok <- !is.na(d$num1) & !is.na(d$ord2)
  expect_equal(R["num1", "ord2"], polyserial(d$num1[ok], d$ord2[ok]))
  expect_equal(R["num1", "num2"], stats::cor(d$num1, d$num2, use = "complete.obs"))
})

test_that("het_cor warns once about unordered factors", {
  d <- mixed_data()
  d$ord1 <- factor(d$ord1, ordered = FALSE)
  expect_warning(
    het_cor(d),
    "Treating unordered factor(s) as ordinal in level order: ord1",
    fixed = TRUE
  )
  expect_no_warning(het_cor(mixed_data()))
})

test_that("het_cor rejects columns that are neither numeric nor factor", {
  d <- mixed_data()
  d$id <- as.character(seq_len(nrow(d)))
  expect_error(het_cor(d), "Columns must be numeric or factor: id")
  expect_error(het_cor(d["num1"]), "at least two columns")
})

test_that("het_cor(plot = TRUE) returns the matrix and a ggplot", {
  res <- het_cor(mixed_data(), plot = TRUE, labels = c(num1 = "First"))
  expect_named(res, c("cor", "plot"))
  expect_equal(res$cor, het_cor(mixed_data()))
  expect_s3_class(res$plot, "ggplot")
})

test_that("plot_cor_triangle shows the lower triangle with type marks", {
  d <- mixed_data()
  p <- plot_cor_triangle(het_cor(d), types = cor_types(d))
  expect_equal(nrow(p$data), 6)
  expect_true(any(grepl("\\*\\*$", p$data$label)))
  expect_true(any(grepl("[0-9]\\*$", p$data$label)))
  expect_match(p$labels$caption, "** Polychoric", fixed = TRUE)
})

test_that("plot_cor_triangle checks labels", {
  R <- het_cor(mixed_data())
  expect_error(plot_cor_triangle(R, labels = c(nope = "x")), "`labels` must be named")
  expect_error(plot_cor_triangle(R, labels = "unnamed"), "`labels` must be named")
})
