
library( testthat)

####testing makeMesh

f <- system.file("ex/test.grd", package="terra")
r <- terra::rast(f)
terra::values( r)[ !is.na( terra::values( r))] <- 1
rm( f)
dat <- simulateData.isdm( expected.pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
terra::crs( dat$covarBrick) <- terra::crs( r)

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
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="SpatRaster")
    testthat::expect_type( object=meshy$n, "integer")
    
    #suppress plotting
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=FALSE)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="SpatRaster")
    testthat::expect_type( object=meshy$n, "integer")
    
    #as I would usually use it
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="SpatRaster")
    testthat::expect_type( object=meshy$n, "integer")
    
    #extensions / finer control
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500, expans.mult=pi/3, max.edge=c(1,0.5), cutoff=c(0.1,0.2), offset=0.91)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="SpatRaster")
    testthat::expect_type( object=meshy$n, "integer")
    
    #extensions / finer control
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500, expandRegion = FALSE)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="SpatRaster")
    testthat::expect_type( object=meshy$n, "integer")
    
    #with some holes
    r <- dat$covarBrick[[1]]
    #values( r) <- 1:ncell( r)
    terra::values(r)[terra::crds(r,na.rm=FALSE)[,1]>179900 & terra::crds(r,na.rm=FALSE)[,1]<180200 & 
                       terra::crds(r,na.rm=FALSE)[,2]>331800 & terra::crds(r,na.rm=FALSE)[,2]<332100] <- NA
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=TRUE, max.n=250, dep.range=500)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="SpatRaster")
    testthat::expect_type( object=meshy$n, "integer")
  }
)

