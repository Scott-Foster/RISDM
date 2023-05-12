
library( testthat)

####testing simulateData.isdm

testthat::test_that(
  "Checking the simulation of data for ISDM examples/checks",
  {
    #using all defaults
    dat <- simulateData.isdm( control=list( doPlot=FALSE))
    testthat::expect_s3_class(object=dat, class="simISDMdata")  #make sure an object has been returned.
    testthat::expect_equal(object=dat$expected.pop.size, expected=4e+5)  #the default
    
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
  }
)

