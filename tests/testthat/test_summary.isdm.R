
library( testthat)

####testing isdm

RandomFields::RFoptions( install="no")
f <- system.file("external/test.grd", package="raster")
r <- raster(f)
raster::values( r)[ !is.na( raster::values( r))] <- 1
rm( f)
dat <- simulateData.isdm( expected.pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
raster::crs( dat$covarBrick) <- raster::crs( r)
meshy <- makeMesh( dat$covarBrick[[1]], max.n=c(500, 150), dep.range=25, expans.mult=20, offset=500, max.edge=5, doPlot=FALSE)
fm <- isdm( observationList=list( POdat=as.data.frame( dat$PO), 
                                        DCdat=as.data.frame( dat$DC),
                                        AAdat=as.data.frame( dat$AA)),
                  covarBrick=dat$covarBrick, 
                  mesh=meshy,
                  responseNames=c( AA="AA"),#, PA="PA"),
                  sampleAreaNames=c( PO=NULL, DC="transectArea", AA="transectArea"),#, PA="transectArea"),
                  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=~1+dist2City,
                  artefactFormulas=list( DC=~1+Survey, AA=~1),#, PA=~1),
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE, 
                                DCmethod="TalyorsLinApprox"))

testthat::test_that(
  "Checking the summary for isdm object (and its print)",
  {
    testthat::expect_s3_class( summary( fm), class="summary.isdm")  #checking summary.isdm
    tmp <- summary( fm)
    testthat::expect_invisible( print( tmp))  #checking print.summary.isdm

    testthat::expect_invisible( print( tmp, digits=6))
    testthat::expect_invisible( print( tmp, digits=0))  #silly
  }
)