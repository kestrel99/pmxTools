test_that("nhanes_peds has the documented structure", {
  expect_s3_class(nhanes_peds, "data.frame")
  expect_named(
    nhanes_peds,
    c("CYCLE", "SEQN", "AGE_MONTHS", "AGE", "SEX", "WT", "HT", "MEC_WT")
  )
  expect_setequal(unique(nhanes_peds$CYCLE), c("2013-14", "2015-16", "2017-18", "2021-23"))
  expect_setequal(unique(nhanes_peds$AGE), 2:17)
  expect_setequal(unique(nhanes_peds$SEX), c("Male", "Female"))
  expect_false(anyNA(nhanes_peds$MEC_WT))
  expect_true(all(nhanes_peds$MEC_WT > 0))
})

test_that("nhanes_peds keeps children with only one of weight or height", {
  has_wt <- !is.na(nhanes_peds$WT)
  has_ht <- !is.na(nhanes_peds$HT)
  expect_true(all(has_wt | has_ht))
  expect_true(any(has_wt & !has_ht))
  expect_true(all(nhanes_peds$WT[has_wt] > 0))
  expect_true(all(nhanes_peds$HT[has_ht] > 0))
})
