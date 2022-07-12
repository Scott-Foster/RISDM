
library( testthat)

####testing isdm

RandomFields::RFoptions( install="no")
f <- system.file("external/test.grd", package="raster")
r <- raster(f)
values( r)[ !is.na( values( r))] <- 1
rm( f)
dat <- simulateData.isdm( expected.pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
crs( dat$covarBrick) <- crs( r)
meshy <- makeMesh( dat$covarBrick[[1]], max.n=c(500, 150), dep.range=25, expans.mult=20, offset=500, max.edge=5, doPlot=FALSE)

testthat::test_that(
  "Checking the model fitting ISDM on irregular boundary",
  {
    fm1 <- list()
    
    #with the PO data only using plugin estimates for now.
    
    fm1[[1]] <- isdm( observationList=list( POdat=as.data.frame( dat$PO)),
                      covarBrick=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( PO="anything"),
                      sampleAreaNames=c( PO=NULL),
                      distributionFormula=~0+Altitude+Temperature,
                      biasFormula=~1+dist2City,
                      artefactFormulas=list( PO=~1),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TalyorsLinApprox"))
    
    fm1[[2]] <- isdm( observationList=list( DCdat=as.data.frame( dat$DC)),
                      covarBrick=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( DC="somebloodything"),
                      sampleAreaNames=c( DC="transectArea"),
                      DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
                      distributionFormula=~0+Altitude+Temperature,
                      artefactFormulas=list( DC=~1+Survey),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TalyorsLinApprox"))
    
    fm1[[3]] <- isdm( observationList=list( AAdat=as.data.frame( dat$AA)),
                      covarBrick=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( AA="AA"),
                      sampleAreaNames=c( AA="transectArea"),
                      distributionFormula=~0+Altitude+Temperature,
                      artefactFormulas=list( AA=~1),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TalyorsLinApprox"))
    
    fm1[[4]] <- isdm( observationList=list( PAdat=as.data.frame( dat$PA)),
                      covarBrick=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( PA="PA"),
                      sampleAreaNames=c( PA="transectArea"),
                      distributionFormula=~0+Altitude+Temperature,
                      artefactFormulas=list( PA=~1),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TalyorsLinApprox"))

    fm1[[5]] <- isdm( observationList=list( POdat=as.data.frame( dat$PO), 
                                            DCdat=as.data.frame( dat$DC),
                                            AAdat=as.data.frame( dat$AA)),
                      #                                        PAdat=as.data.frame( simDat1$PA)),
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
  }
)
