library(testthat)

test_that("gcv_convert", {
  expect_equal(
    gcv_convert(0.01),
    100*sqrt(exp(0.01)-1),
    info="The math is numerically accurate"
  )
  expect_equal(
    gcv_convert(c(0.1, 0.01)),
    100*sqrt(exp(c(0.1, 0.01))-1),
    info="Vectorized input works"
  )
  expect_equal(
    gcv_convert(0.01),
    gcv_convert(gsd=0.1),
    info="gvar and gsd work the same"
  )
  expect_error(
    gcv_convert(gvar=0.01, gsd=0.1),
    regexp="Only one of `gvar` or `gsd` may be provided at a time.",
    fixed=TRUE
  )
})
