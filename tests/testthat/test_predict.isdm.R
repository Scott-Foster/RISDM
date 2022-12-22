
library( testthat)

####testing isdm

f <- system.file("external/test.grd", package="raster")
r <- raster::raster(f)
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
  "Checking the prediction from an isdm object.  Predicting to original raster only.",
  {
    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=500)
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")
    tmp <- raster::stack( raster::crop( dat$covarBrick$Intensity, fm$preds$mean.field), fm$preds$mean.field)
    raster::plot( tmp, nc=2)

    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=50, n.batches = 3)
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")
    
    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5)  #checking another value of S
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")

    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5, intercept.terms="Intercept.DC:Survey_DCDCsurvey_2")
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")

    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5, n.threads=4)
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")

    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5, includeRandom=FALSE)
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")

    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5, includeFixed=FALSE)
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")

    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5, includeBias=FALSE)
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")
    
    fm$preds <- predict( object=fm, covarRaster=dat$covarBrick, S=5, type="probability")
    testthat::expect_s4_class( fm$preds$mean.field, class="RasterBrick")
  }
)
