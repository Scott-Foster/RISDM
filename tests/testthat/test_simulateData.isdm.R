
library( testthat)

####testing simulateData.isdm

testthat::test_that(
  "Checking the simulation of data for ISDM examples/checks",
  {
    #using all defaults
    dat <- simulateData.isdm( control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=dat$expected.pop.size, expected=10000)  #the default
    
    #setting expected.pop.size
    dat <- simulateData.isdm( expected.pop.size = 5000, control=list(doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=dat$expected.pop.size, expected=5000)  #the default
    
    #setting expectred.n.PO
    dat <- simulateData.isdm( expected.n.PO=5, control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
#    testthat::expect_equal(object=nrow( sp::coordinates( dat$PO)), expected=5, tolerance=0.75)  #wide tolerance (abs( x-y) / abs( y)) < tolerance
    
    #setting n.PA, n.AA and n.DC
    dat <- simulateData.isdm( n.PA=25, n.AA=1053, n.DC=74, control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=as.numeric( sapply( dat[c("PA","AA","DC")], nrow)), expected=c(25,1053,74))
    
    #a raster to test
    r <- terra::rast( system.file( "extdata/ACT_DemoData.grd", package="RISDM"))
    r <- terra::aggregate( r, 3) #to make it finite time (if not using fft approach)
    n1 <- sum( !is.na( terra::values( r[[1]])))
    #giving raster boundary
    dat <- simulateData.isdm( rasterBoundary=r[[1]], control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=sum( !is.na( terra::values( dat$covarBrick$biasIntensity))), expected=n1)
    #giving our own covariates and bias covariates
    dat <- simulateData.isdm( rasterCovars=r[[2:3]], rasterBiasCovar = r[[1]], control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=sum( !is.na( terra::values( dat$covarBrick$biasIntensity))), expected=n1)
  }
)

