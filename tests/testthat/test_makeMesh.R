
library( testthat)

####testing makeMesh

f <- system.file("external/test.grd", package="raster")
r <- raster::raster(f)
raster::values( r)[ !is.na( raster::values( r))] <- 1
rm( f)
dat <- simulateData.isdm( expected.pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
raster::crs( dat$covarBrick) <- raster::crs( r)

#dat <- simulateData.isdm( control=list(doPlot=FALSE))

testthat::test_that(
  "checking creation of the mesh",
  {
    #this shouldn't be needed, but the INLA code seems to spark a warning from Matrix (only first run). 
    #Let's get the first run out of the way...
    meshy <- makeMesh( dat$covarBrick[[1]])
    
    #using all defaults
    meshy <- makeMesh( dat$covarBrick[[1]])
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
    
    #suppress plotting
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=FALSE)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
    
    #as I would usually use it
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
    
    #extensions / finer control
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500, expans.mult=pi/3, max.edge=c(1,0.5), cutoff=c(0.1,0.2), offset=0.91)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
    
    #extensions / finer control
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500, expandRegion = FALSE)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
  }
)

