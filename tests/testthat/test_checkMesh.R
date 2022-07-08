
library( testthat)

####testing makeMesh

RandomFields::RFoptions( install="no")
dat <- simulateData.isdm( control=list(doPlot=FALSE))
meshy <- makeMesh( dat$covarBrick[[1]])

testthat::test_that(
  "testing code for checking meshes",
  {
    testthat::expect_invisible(checkMesh( meshy, meshy$hull))
  }
)

