
library( testthat)

####testing residual plotting for isdm

f <- system.file("ex/test.grd", package="terra")
r <- terra::rast(f)
terra::values( r)[ !is.na( terra::values( r, na.rm=FALSE))] <- 1
rm( f)
dat <- simulateData.isdm( expected.pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
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
            distributionFormula=~0+Altitude+Temperature,
            biasFormula=~1+dist2City,
            artefactFormulas=list( DC=~1+Survey, AA=~1),#, PA=~1),
            control=list( int.prec=0.01, other.prec=1,
                          calcICs=FALSE,
                          prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                          coord.names=c("x","y"),
                          n.threads=8,
                          addRandom=TRUE, 
                          DCmethod="TaylorsLinApprox"))

testthat::test_that(
  "Checking the production of residuals (for diagnostics) for isdm object.",
  {
    my.resids <- residuals( fm, covars=dat$covarBrick)
    testthat::expect_type( object=my.resids, type="list")
  }
)


