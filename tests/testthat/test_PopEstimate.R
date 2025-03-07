
library( testthat)

####testing isdm

f <- system.file("ex/test.grd", package="terra")
r <- terra::rast(f)
terra::values( r)[ !is.na( terra::values( r, na.rm=FALSE))] <- 1
rm( f)
dat <- simulateData.isdm( pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
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
            artefactFormulas=list( DC=~1, AA=~1),#, PA=~1),
            control=list( int.prec=0.01, other.prec=1,
                          calcICs=FALSE,
                          prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
                          coord.names=c("x","y"),
                          n.threads=8,
                          addRandom=TRUE, 
                          DCmethod="TaylorsLinApprox"))

fm$preds <- predict( object=fm, covars=dat$covarBrick, S=50)

testthat::test_that(
  "Checking the prediction from an isdm object.  Predicting to original raster only.",
  {
    fm$popEsti <- PopEstimate( preds=fm$preds)
    testthat::expect_length( fm$popEsti, 6)
    
    tmp <- PopEstimate( preds=fm$preds, control=list( probs=c(0.0025,0.4,0.6,0.9975)))
    testthat::expect_length( tmp, 6)
    
    tmp <- PopEstimate( preds=fm$preds, intercept.terms = c( "DC_Intercept","DC_SurveyDCsurvey_3"))
    testthat::expect_length( tmp, 6)
  }
)
    
