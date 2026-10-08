test_that("weighted_quantile is the type 4 quantile for equal weights", {
  x <- c(5, 1, 3, 2, 4, 8, 6)
  p <- c(0.1, 0.25, 0.5, 0.9)
  expect_equal(
    weighted_quantile(x, rep(1, 7), p),
    stats::quantile(x, p, type = 4, names = FALSE)
  )
})

test_that("weighted_quantile follows the weights", {
  # Cumulative weights 0.25, 0.5, 1 at x = 1, 2, 3, interpolated linearly
  x <- c(1, 2, 3)
  w <- c(1, 1, 2)
  expect_equal(weighted_quantile(x, w, c(0.5, 0.75, 1)), c(2, 2.5, 3))
  expect_equal(weighted_quantile(x, w * 10, 0.75), 2.5)
})

test_that("weighted_sd is the population SD for equal weights", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  expect_equal(weighted_sd(x, rep(1, 8)), 2)
})

test_that("weighted_cor matches cor() for equal weights", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 1, 4, 3, 5)
  expect_equal(weighted_cor(x, y, rep(1, 5)), cor(x, y))
})

test_that("robust_scale is min(SD, IQR/1.349) and falls back to 0", {
  x <- c(1, 2, 3, 4, 100)
  w <- rep(1, 5)
  iqr <- diff(weighted_quantile(x, w, c(0.25, 0.75))) / 1.349
  expect_equal(robust_scale(x, w), min(weighted_sd(x, w), iqr))
  expect_equal(robust_scale(c(3, 3, 3), rep(1, 3)), 0)
})
