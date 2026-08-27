library(vdiffr)

set.seed(740727)

df <- data.frame(wt = c(rnorm(1000, mean=70, sd=15), rnorm(1000,mean=85, sd=20)), 
                 Group=rep(c("Group 1","Group 2"), each=1000))

# SVG snapshots record text metrics from the system "sans" font, which resolves
# differently on macOS (Helvetica) and Windows (Arial) than on the Linux machine
# the snapshots are generated on. Compare figures on Linux only.
test_that("plot_dist", {
  skip_on_os(c("mac", "windows"))
  expect_doppelganger("Single distplot", plot_dist(df[df$Group=="Group 1",], yvar="wt", ylb="Weight (kg)"))
  expect_doppelganger("Conditioned distplot", plot_dist(df, yvar="wt", xvar="Group", ylb="Weight (kg)"))
})

test_that("plot_nmprogress", {
  skip_on_os(c("mac", "windows"))
  expect_doppelganger("Est", plot_nmprogress(fileName = "testData/runR001", fileExt = ".res", metric = "est"))
  expect_doppelganger("Perc", plot_nmprogress(fileName = "testData/runR001", fileExt = ".res"))
})