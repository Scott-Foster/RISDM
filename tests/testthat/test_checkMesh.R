
library( testthat)

####testing makeMesh

f <- system.file("external/test.grd", package="raster")
r <- raster::raster(f)

dat <- simulateData.isdm( control=list(doPlot=FALSE))
raster::crs( dat$covarBrick) <- raster::crs(r)
meshy <- makeMesh( dat$covarBrick[[1]])

raster::values( r)[ !is.na( raster::values( r))] <- 1
rm( f)
dat1 <- simulateData.isdm( expected.pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
raster::crs( dat1$covarBrick) <- raster::crs( r)
meshy1 <- makeMesh( dat1$covarBrick[[1]], max.n=c(500, 150), dep.range=25, expans.mult=20, offset=500, max.edge=5, doPlot=FALSE)

testthat::test_that(
  "testing code for checking meshes",
  {
    #avoids the warning thrown by not having a CRS for the data.
    testthat::expect_warning(checkMesh( meshy, meshy$hull))
    #shouldn't produce an error as the CRS is defined
    testthat::expect_invisible( checkMesh( meshy1, meshy1$hull))
  }
)

