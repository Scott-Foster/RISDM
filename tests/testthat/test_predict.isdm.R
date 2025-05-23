
library( testthat)

####testing isdm

f <- system.file("ex/test.grd", package="terra")
r <- terra::rast(f)
terra::values( r)[ !is.na( terra::values( r, na.rm=FALSE))] <- 1
rm( f)
set.seed( 747)
dat <- simulateData.isdm( pop.size=5000, transect.size=0.6, rasterBoundary=r, control=list(doPlot=FALSE))
terra::crs( dat$covarBrick) <- terra::crs( r)
meshy <- makeMesh( dat$covarBrick[[1]], max.n=c(500, 150), dep.range=25, expans.mult=20, offset=500, max.edge=5, doPlot=FALSE)

fm <- isdm( observationList=list( POdat=as.data.frame( dat$PO), 
                                  DCdat=as.data.frame( dat$DC),
                                  AAdat=as.data.frame( dat$AA)),
            covars=dat$covarBrick, 
            mesh=meshy,
            responseNames=c( AA="AA"),#, PA="PA"),
            sampleAreaNames=c( PO=NULL, DC="transectArea", AA="transectArea"),#, PA="transectArea"),
            DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
            distributionFormula=~0+var1,
            biasFormula=~1+biasLinPred,
            artefactFormulas=list( DC=~1+Survey, AA=~1),#, PA=~1),
            control=list( int.prec=0.01, other.prec=1,
                          calcICs=FALSE,
                          prior.range=c( 25,0.1), prior.space.sigma=c( 2.5,0.1),
                          coord.names=c("x","y"),
                          n.threads=8,
                          addRandom=TRUE, 
                          DCmethod="TaylorsLinApprox"))

fm1 <- isdm( observationList=list( POdat=as.data.frame( dat$PO), 
                                  DCdat=as.data.frame( dat$DC),
                                  AAdat=as.data.frame( dat$AA)),
            covars=dat$covarBrick, 
            mesh=meshy,
            responseNames=c( AA="AA"),#, PA="PA"),
            sampleAreaNames=c( PO=NULL, DC="transectArea", AA="transectArea"),#, PA="transectArea"),
            DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
            distributionFormula=~0+var1,
            biasFormula=~1+biasLinPred,
            artefactFormulas=list( DC=~1+Survey, AA=~1),#, PA=~1),
            control=list( int.prec=0.01, other.prec=1,
                          calcICs=FALSE,
                          prior.range=c( 25,0.1), prior.space.sigma=c( 2.5,0.1),
                          coord.names=c("x","y"),
                          n.threads=8,
                          addRandom=FALSE, 
                          DCmethod="TaylorsLinApprox"))


testthat::test_that(
  "Checking the prediction from an isdm object.  Predicting to original raster only.",
  {
    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=500)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
    tmp <- c( terra::crop( dat$covarBrick$Intensity, fm$preds$field), fm$preds$field)
    terra::plot( tmp, nc=2)

    #checking the confidence levels
    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=500, confidence.level=0.8)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
    tmp <- c( terra::crop( dat$covarBrick$Intensity, fm$preds$field), fm$preds$field)
    terra::plot( tmp, nc=2)
    
    fm1$preds <- predict( object=fm1, covars=dat$covarBrick, S=500)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
    
    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=50, n.batches = 3)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
    
    fm1$preds <- predict( object=fm1, covars=dat$covarBrick, S=500, n.batches = 3)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
    
    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5)  #checking another value of S
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")

    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5, intercept.terms=c( "DC_Intercept", "DC_Survey"))
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")

    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5, n.threads=4)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")

    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5, includeRandom=FALSE)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")

    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5, includeFixed=FALSE)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")

    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5, includeBias=FALSE)
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
    
    fm$preds <- predict( object=fm, covars=dat$covarBrick, S=5, type="probability")
    testthat::expect_s4_class( fm$preds$field, class="SpatRaster")
  }
)
