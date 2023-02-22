library(vdiffr)

set.seed(740727)

df <- data.frame(wt = c(rnorm(1000, mean=70, sd=15), rnorm(1000,mean=85, sd=20)), 
                 Group=rep(c("Group 1","Group 2"), each=1000))

test_that("plot_dist", {
  expect_doppelganger("Single distplot", plot_dist(df[df$Group=="Group 1",], yvar="wt", ylb="Weight (kg)"))
  expect_doppelganger("Conditioned distplot", plot_dist(df, yvar="wt", xvar="Group", ylb="Weight (kg)"))
})

test_that("plot_nmprogress", {
  expect_doppelganger("Est", plot_nmprogress(fileName = "testData/runR001", fileExt = ".res", metric = "est"))
  expect_doppelganger("Perc", plot_nmprogress(fileName = "testData/runR001", fileExt = ".res"))
})