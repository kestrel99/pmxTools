test_that("read_nm", {
  expect_equal(
    read_nm("testData/runR001.xml", quiet = TRUE),
    readRDS("testData/read_nm.RDS")
  )
})

test_that("rnm", {
  expect_equal(
    rnm(index = "001", prefix="runR", pathNM = "testData", ext=".res", extmod = ".ctl"),
    readRDS("testData/rnm.RDS")
  )
})

test_that("read_nmcov", {
  expect_equal(
    read_nmcov("testData/runR001.cov", quiet = TRUE),
    readRDS("testData/nmcov.RDS")
  )
})

test_that("read_nmext", {
  expect_equal(
    read_nmext("testData/runR001.ext", quiet = TRUE),
    readRDS("testData/nmext.RDS")
  )
})

nm1 <- read_nm("testData/runR001.xml", quiet = TRUE)

test_that("get_theta", {
  expect_equal(
    get_theta(nm1),
    c(THETA1=3.980130, THETA2=68.219400, THETA3=0.199472),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_theta(nm1, output="se"),
    c(THETA1=0.09897590, THETA2=1.92949000, THETA3=0.00280475),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_theta(nm1, output="rse"),
    c(THETA1=2.48675, THETA2=2.82836, THETA3=1.40609),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_theta(nm1, output="95ci", sigdig = 3),
    c(THETA1="3.79-4.17", THETA2="64.4-72.0", THETA3="0.194-0.205"),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_theta(nm1, output="all", sigdig = 4),
    list(
      Theta = c(THETA1=3.9800, THETA2=68.2200, THETA3=0.1995),
      ThetaSE = c(THETA1=0.098980, THETA2=1.929000, THETA3=0.002805),
      ThetaRSE = c(THETA1=2.487, THETA2=2.828, THETA3=1.406),
      Theta95CI = c(THETA1="3.786-4.174", THETA2="64.44-72.00", THETA3="0.1940-0.2050")
    ),
    info="The function returns the expected answer"
  )
})

test_that("get_omega", {
  expect_equal(
    get_omega(nm1),
    matrix(data=c(0.0715555,0.0000000,0.0000000,0.0921585), nrow = 2, ncol = 2,
           dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_omega(nm1, output="se"),
    matrix(data=c(0.00912817,NA,NA,0.0122009), nrow = 2, ncol = 2,
           dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_omega(nm1, output="rse"),
    matrix(data=c(12.7568,NA,NA,13.239), nrow = 2, ncol = 2,
           dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_omega(nm1, output="cor"),
    matrix(data=c(0.267499,0,0,0.303576), nrow = 2, ncol = 2,
           dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_omega(nm1, output="cse"),
    matrix(data=c(0.0170621,NA,NA,0.0200952), nrow = 2, ncol = 2,
           dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_omega(nm1, output="95ci", sigdig = 3),
    matrix(data=c("0.0537-0.0894",NA,NA,"0.0682-0.116"), nrow = 2, ncol = 2,
           dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_omega(nm1, output="all", sigdig = 4),
    list(
      Omega = matrix(data=c(0.07156,0.0000000,0.0000000,0.09216), nrow = 2, ncol = 2,
                     dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
      OmegaSE = matrix(data=c(0.009128,NA,NA,0.0122), nrow = 2, ncol = 2,
                       dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
      OmegaRSE = matrix(data=c(12.76,NA,NA,13.24), nrow = 2, ncol = 2,
                        dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
      Omega95CI = matrix(data=c("0.05366-0.08945",NA,NA,"0.06824-0.1161"), nrow = 2, ncol = 2,
                         dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
      OmegaCorrelation = matrix(data=c(0.2675,0,0,0.3036), nrow = 2, ncol = 2,
                       dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2"))),
      OmegaCorrelationSE = matrix(data=c(0.01706,NA,NA,0.0201), nrow = 2, ncol = 2,
                        dimnames = list(c("OMEGA1","OMEGA2"), c("OMEGA1","OMEGA2")))
    ),
    info="The function returns the expected answer"
  )
})

test_that("get_shrinkage", {
  expect_equal(
    get_shrinkage(nm1),
    c(ETA1=1.37, ETA2=1.46),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_shrinkage(nm1,output = "epsilon"),
    c(EPS1=5.22),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_shrinkage(nm1,output = "all"),
    list(
      etasd = c(ETA1=1.37, ETA2=1.46),
      etavr = c(ETA1=2.73, ETA2=2.90),
      epsilonsd = c(EPS1=5.22),
      epsilonvr = c(EPS1=10.2)
    ),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_shrinkage(nm1, type = "vr", sigdig = 4),
    c(ETA1=2.725, ETA2=2.902),
    info="The function returns the expected answer"
  )
  expect_equal(
    get_shrinkage(nm1,output = "epsilon", type = "vr", sigdig = 4),
    c(EPS1=10.17),
    info="The function returns the expected answer"
  )

})

test_that("get_est_table", {
  expect_equal(
    get_est_table(nm1, thetaLabels = c("CL","V","KA"), omegaLabels = c("CL","V"), sigmaLabels = "Residual"),
    readRDS("testData/get_est_table.RDS")
  )
})

test_that("get_probinfo", {
  expect_equal(
    get_probinfo(nm1),
    readRDS("testData/get_probinfo.RDS")
  )
})

test_that("sample_uncertainty", {
  expect_equal(
    suppressMessages(sample_uncert("testData/runR001.xml", n=50, seed=740727)),
    readRDS("testData/sample_uncert.RDS")
  )
})

test_that("sample_omega", {
  expect_equal(
    suppressMessages(sample_omega("testData/runR001.xml", n=50, seed=740727)),
    readRDS("testData/sample_omega.RDS")
  )
})

test_that("sample_sigma", {
  expect_equal(
    suppressMessages(sample_sigma("testData/runR001.xml", n=50, seed=740727)),
    readRDS("testData/sample_sigma.RDS")
  )
})
