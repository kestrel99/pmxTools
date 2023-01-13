context("blq_trans")

test_that("estimate_lloq works", {
  expect_equal(estimate_lloq(2), 2)
  expect_equal(estimate_lloq(c(NA, 2)), 2)
  expect_equal(estimate_lloq(c(NA, -1, 2)), 2)
  expect_equal(
    expect_warning(
      estimate_lloq(c(NA, -1)),
      regexp="No samples above the lloq, using 1",
      fixed=TRUE
    ),
    1
  )
})

test_that("blq_trans forward and inverse", {
  expect_equal(
    ftrans_blq_linear(lloq=1, multiplier=0.5)(c(0, 2)),
    c(0.5, 2)
  )
  expect_equal(
    ftrans_blq_linear(lloq=2, multiplier=0.5)(c(0, 2)),
    c(1, 2)
  )
  expect_equal(
    ftrans_blq_linear(lloq=3, multiplier=0.8)(c(0, 2)),
    c(2.4, 2.4)
  )

  expect_equal(
    ftrans_blq_log(lloq=1, multiplier=0.5, base=10)(c(0, 2)),
    log10(c(0.5, 2))
  )
  expect_equal(
    ftrans_blq_log(lloq=2, multiplier=0.5)(c(0, 2)),
    log10(c(1, 2))
  )
  expect_equal(
    ftrans_blq_log(lloq=3, multiplier=0.8)(c(0, 2)),
    log10(c(2.4, 2.4))
  )
  expect_equal(
    ftrans_blq_log(lloq=3, multiplier=0.8, base=2)(c(0, 2)),
    log(c(2.4, 2.4), base=2)
  )
  
  expect_equal(
    itrans_blq_linear(0.1)(c(0, 2)),
    c(0.1, 2)
  )
  expect_equal(
    itrans_blq_log(0.1, base=10)(c(-10, 0, 2)),
    c(0.1, 1, 100)
  )
})

test_that("blq breaks", {
  expect_equal(
    breaks_blq_general(lloq=3, breakfun=scales::breaks_extended)(NA_real_, n=5),
    scales::breaks_extended(n=5)(NA_real_),
    info="NA falls back to the normal method (whatever that may be)"
  )
  expect_equal(
    breaks_blq_general(lloq=3, breakfun=scales::breaks_extended)(c(1, NA), n=5),
    3,
    info="All BLQ or NA returns just the lloq"
  )
  expect_equal(
    breaks_blq_general(lloq=3, breakfun=scales::breaks_extended)(c(1:3, 3.0001), n=5),
    scales::breaks_extended(n=5)(c(3, 3.0001)),
    info="If data are only slightly above the LLOQ, then the scales are made from the LLOQ to the top of the data."
  )
  expect_equal(
    breaks_blq_general(lloq=3, breakfun=scales::breaks_extended)(1:100, n=5),
    c(3, 25, 50, 75, 100)
  )
  expect_equal(
    breaks_blq_general(lloq=3, breakfun=scales::breaks_log, trans=log)(1:100, n=5),
    c(3, 10, 30, 100, 300),
    info="The scale is created only on data above the LLOQ (note that the result is different with just scales::breaks_log()(1:100)."
  )
  expect_equal(
    breaks_blq_general(lloq=9, breakfun=scales::breaks_log, trans=log)(1:100, n=5),
    c(9, 30, 100, 300),
    info="Values that are too close together are dropped"
  )
  expect_equal(
    breaks_blq_general(lloq=9, breakfun=scales::breaks_log, trans=log)(c(1, 90:100), n=5),
    c(9, 30, 100, 300),
    info="A single value BLQ will trigger the lloq to be in the scale"
  )
  expect_equal(
    breaks_blq_general(lloq=9, breakfun=scales::breaks_log, trans=log)(90:100, n=5),
    scales::breaks_log(n=5)(90:100),
    info="No data BLQ means that the lloq will not be on the scale"
  )
})

test_that("label_blq works", {
  expect_equal(
    label_blq(lloq=2)(1:3),
    c("<2", "<2", "3"),
    info="lloq_text is imputed as the <LLOQ"
  )
  expect_equal(
    label_blq(lloq=2, lloq_text="BLQ")(1:3),
    c("BLQ", "BLQ", "3"),
    info="lloq_text is imputed as the <LLOQ"
  )
  expect_equal(
    label_blq(lloq=2, lloq_text="BLQ")(2*(1+sqrt(.Machine$double.eps)/2)),
    "BLQ",
    info="lloq_text is used, even if very nearly LLOQ (to allow for machine precision differences)"
  )
})
