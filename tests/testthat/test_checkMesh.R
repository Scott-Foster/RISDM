
library( testthat)

####testing makeMesh

f <- system.file("ex/test.grd", package="terra")
r <- terra::rast(f)
rm( f)

meshy1 <- makeMesh( r, max.n=c(1000, 250), dep.range=25, expans.mult=20, offset=500, max.edge=5, doPlot=TRUE)

testthat::test_that(
  "testing code for checking meshes",
  {
    testthat::expect_invisible( checkMesh( meshy1))
  }
)

