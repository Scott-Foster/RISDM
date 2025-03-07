
library( testthat)

####testing simulateData.isdm

testthat::test_that(
  "Checking the simulation of data for ISDM examples/checks",
  {
    #using all defaults
    dat <- simulateData.isdm( control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=dat$pop.size, expected=10000)  #the default
    
    #setting expected.pop.size
    dat <- simulateData.isdm( pop.size = 5000, control=list(doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=dat$pop.size, expected=5000)  #the default
    
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
    dat <- simulateData.isdm( pop.size=10000, distForm= ~soilMoisture+MeanMinTemp, biasForm= ~1+logAccessibility, covarBrick=r, control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=sum( !is.na( terra::values( dat$covarBrick$Intensity))), expected=n1)
    
    #giving an intensity
    dat <- simulateData.isdm( pop.size=10000, biasForm=~1+logAccessibility, covarBrick=r, Intensity=exp( r$soilMoisture), control=list(doPlot=FALSE))
    testthat::expect_s3_class( object=dat, class="simISDMdata")
    
    #giving own locations for surveys
    tmp <- terra::as.data.frame( r, xy=TRUE, cells=TRUE, na.rm=TRUE)
    tmp.pa <- as.matrix( tmp[sample.int(n=nrow(tmp), size=50, replace=FALSE),c("x","y")])
    tmp.aa <- as.matrix( tmp[sample.int(n=nrow(tmp), size=25, replace=FALSE),c("x","y")])
    tmp.dc <- as.matrix( tmp[sample.int(n=nrow(tmp), size=30, replace=FALSE),c("x","y")])
    
    dat <- simulateData.isdm( covarBrick = r, biasForm=~1+logAccessibility, Intensity=exp(r$soilMoisture),
                       PA.locats=tmp.pa[,c("x","y")], AA.locats=tmp.aa[,c("x","y")], DC.locats=tmp.dc[,c("x","y")], control=list( doPlot=FALSE))
    testthat::expect_s3_class( object=dat, class="simISDMdata")
    testthat::expect_equal(object=as.numeric( sapply( dat[c("PA","AA","DC")], nrow)), expected=c(nrow( tmp.pa),nrow(tmp.aa),nrow(tmp.dc)))

  }
)

