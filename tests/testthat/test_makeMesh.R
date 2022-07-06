
library( testthat)

####testing makeMesh

RandomFields::RFoptions( install="no")
dat <- simulateData.isdm( control=list(doPlot=FALSE))

testthat::test_that(
  "checking creation of the mesh",
  {
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
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=FALSE, max.n=250, dep.range=2)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
    
    #extensions / finer control
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=FALSE, max.n=250, dep.range=2, expandRegion = FALSE)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
    
    meshy <- makeMesh( dat$covarBrick[[1]], doPlot=FALSE, max.n=250, dep.range=2, expans.mult=pi/3, max.edge=c(1,0.5), cutoff=c(0.1,0.2), offset=0.91)
    testthat::expect_s3_class(object=meshy, class="inla.mesh")  #make sure an object has been returned.
    testthat::expect_s4_class( object=meshy$risdmBoundary$ras$lower.res, class="RasterLayer")
    testthat::expect_type( object=meshy$n, "integer")
  }
)

