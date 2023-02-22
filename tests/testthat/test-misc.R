library(testthat)
library(pmxTools)

### AUC

test_that("AUC", {
  d <- expand.grid(ID = 1:5, TIME=seq(0,24,by=0.5))
  d <- d[order(d$ID, d$TIME),]
  d$DV <- c(calc_sd_1cmt_linear_bolus(t=seq(0,24,by=0.5), CL=6, V=25, dose=600),
            calc_sd_1cmt_linear_bolus(t=seq(0,24,by=0.5), CL=5.2, V=26, dose=600),
            calc_sd_1cmt_linear_bolus(t=seq(0,24,by=0.5), CL=3, V=27, dose=600),
            calc_sd_1cmt_linear_bolus(t=seq(0,24,by=0.5), CL=4.2, V=30, dose=600),
            calc_sd_1cmt_linear_bolus(t=seq(0,24,by=0.5), CL=9, V=20, dose=600))
  a <- get_auc(d)
  a$AUC <- signif(a$AUC,5)

  expect_equal(a, data.frame(ID=1:5, AUC=c(99.804, 114.530, 186.150, 137.950,  66.946)))
})

test_that("count_na", {
  expect_equal(
    count_na(c(0,5,7,NA,3,3,NA)),
    2,
    info="The function returns the expected answer"
  )
  expect_warning(out <- count_na(c(0,5,7,NA,3,3,NA,NaN,Inf,NaN)), "2 NaN values included in the NA count.")
  expect_equal(out, 4)
  
})

test_that("table_rtf", {
  expect_snapshot(table_rtf(iris))
  
})
