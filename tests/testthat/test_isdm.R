
library( testthat)

####testing isdm

f <- system.file("ex/meuse.tif", package="terra")
r <- terra::rast(f)
terra::values( r)[ !is.na( terra::values( r))] <- 1
rm( f)
dat <- simulateData.isdm( pop.size=200000, rasterBoundary=r, control=list(doPlot=FALSE))
terra::crs( dat$covarBrick) <- terra::crs( r)
meshy <- makeMesh( dat$covarBrick[[1]], max.n=c(500, 250), dep.range=25, expans.mult=15, offset=250, max.edge=5, doPlot=FALSE)
#terra::values( r)[terra::crds(r)[,1] > 179900 & terra::crds(r)[,1] < 180100 & terra::crds(r)[,2]>331400 & terra::crds(r)[,2]<331600] <- NA

testthat::test_that(
  "Checking the model fitting ISDM on irregular boundary",
  {
    fm1 <- list()
    
    fm1[[1]] <- isdm( observationList=list( POdat=as.data.frame( dat$PO)),
            covars=dat$covarBrick, 
            mesh=meshy,
            responseNames=NULL,
            sampleAreaNames=NULL,
            distributionFormula=~0+var1,
            biasFormula=~1+biasLinPred,
            artefactFormulas=list( PO=~1),
            control=list( coord.names=c("x","y"), 
                          int.sd=1000, other.sd=10, prior.mean=0,
                          prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
			  n.threads=8, addRandom=TRUE,
			  DCmethod="TaylorsLinApprox"))
			  
    fm1[[2]] <- isdm( observationList=list( DCdat=dat$DC),
                      covars=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( DC="somebloodything"),
                      sampleAreaNames=c( DC="transectArea"),
                      DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
		      distributionFormula=~0+var1,
                      artefactFormulas=list( DC=~1),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TaylorsLinApprox"))
    
    fm1[[3]] <- isdm( observationList=list( AAdat=as.data.frame( dat$AA)),
                      covars=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( AA="AA"),
                      sampleAreaNames=c( AA="transectArea"),
                      distributionFormula=~0+var1,
                      artefactFormulas=list( AA=~1),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TaylorsLinApprox"))
    
    fm1[[4]] <- isdm( observationList=list( PAdat=as.data.frame( dat$PA)),
                      covars=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( PA="PA"),
                      sampleAreaNames=c( PA="transectArea"),
                      distributionFormula=~0+var1,
                      artefactFormulas=list( PA=~1),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TaylorsLinApprox"))

    fm1[[5]] <- isdm( observationList=list( POdat=as.data.frame( dat$PO), 
                                            DCdat=as.data.frame( dat$DC),
                                            AAdat=as.data.frame( dat$AA)),#,
#                                            PAdat=as.data.frame( dat$PA)),
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
                                    prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TaylorsLinApprox"))
				    
    fm1[[6]] <- isdm( observationList=list( DCdat=as.data.frame( dat$DC)),
                      covars=dat$covarBrick, 
                      mesh=meshy,
                      responseNames=c( DC="somebloodything"),
                      sampleAreaNames=c( DC="transectArea"),
                      DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
                      distributionFormula=~0+stats::poly( var1, degree=2),
                      artefactFormulas=list( DC=~1+Survey:stats::poly( var1, degree=2)),
                      control=list( int.prec=0.01, other.prec=1,
                                    calcICs=FALSE,
                                    prior.range=c(25,0.1), prior.space.sigma=c( 2.5,0.1),
                                    coord.names=c("x","y"),
                                    n.threads=8,
                                    addRandom=TRUE, 
                                    DCmethod="TaylorsLinApprox"))
				    
    testthat::expect_length( fm1, 6)
    testthat::expect_s3_class(object=fm1[[1]], class="isdm")
    testthat::expect_s3_class(object=fm1[[2]], class="isdm")
    testthat::expect_s3_class(object=fm1[[3]], class="isdm")
    testthat::expect_s3_class(object=fm1[[4]], class="isdm")
    testthat::expect_s3_class(object=fm1[[5]], class="isdm")
    testthat::expect_s3_class(object=fm1[[6]], class="isdm")
  }
)
