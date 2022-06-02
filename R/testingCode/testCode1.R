
source( "~/ABARES/estimationCode/R/simulate1.R")

my.coefs <- list(dist=c(-1,-0.25,0.75), bias=c(-2,-1.5))

##########################################
####  Is the simulation process OK??

#Abundance (AA) checking
simDat <- simulate.isdm( n.AA=500, coefs=my.coefs)
predDat <- as.data.frame( simDat$covarBrick)
predDat$transectArea <- prod( res( simDat$covarBrick))

tmpDat <- as.data.frame( simDat$AA)
tmpDat <- cbind( tmpDat, extract( simDat$covarBrick, tmpDat[,c("x","y")]))
fm.AA <- glm( AA~1+Altitude+Temperature, family=poisson(), offset=log(transectArea), data=tmpDat)
summary( fm.AA)  #seems to be pretty good!  Upto spatial RE confounding?
pred.AA <- predict( fm.AA, newdata=predDat, type='response')
print( sum( pred.AA) / simDat$expected.pop.size)

#Presence-absence (PA) checking
simDat <- simulate.isdm( n.PA=500, coefs=my.coefs)

tmpDat <- as.data.frame( simDat$PA)
tmpDat <- cbind( tmpDat, extract( simDat$covarBrick, tmpDat[,c("x","y")]))
fm.PA <- glm( PA~1+Altitude+Temperature, family=binomial( link="cloglog"), offset=log( transectArea), data=tmpDat)
summary( fm.PA)  #seems to be pretty good!  Upto spatial RE confounding?
pred.PA <- exp( predict( fm.PA, newdata=predDat, type='link'))
print( sum( pred.AA) / simDat$expected.pop.size)

####  Simulation seems to be appropriate.  Although this doesn't check the PO data

############################################
####  Checking the current github version of RISDM

library( RISDM)
library( raster)
source( "~/ABARES/estimationCode/R/simulate1.R")

simDat <- simulate.isdm( expected.n.PO=3000, n.PA=500, n.AA=500, coefs=my.coefs)

dev.off()
meshy <- RISDM::makeMesh( ras=simDat$covarBrick[[1]], max.n=c(2000,350), dep.range=1)#, expans.mult=2.5)

fm.isdm <- RISDM::isdm( POdat=as.data.frame( simDat$PO), AAdat=as.data.frame( simDat$AA), PA=as.data.frame( simDat$PA),
            covarBrick=simDat$covarBrick, mesh=meshy, boundary=meshy$risdmBoundary$poly,
            responseNames=c(PO="anything", AA="AA", PA="PA"), 
            sampleAreaNames=c( AA="transectArea", PA="transectArea"),
            distributionFormula=~0+Altitude+Temperature,
            biasFormulas=list( PO=~1+dist2City, AA=~1, PA=~1),
            #DCobserverInfo=list( SurveyID="SurveyDC", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
            control=list( int.prec=0.01, other.prec=1,
                          calcICs=FALSE,
                          prior.range=c(0.05,0.1), prior.space.sigma=c( 2,0.0001),
                          coord.names=c("x","y"),
                          n.threads=8))
summary( fm.isdm$mod)
plot( fm.isdm$preds)

##using new code
source( "~/ABARES/estimationCode/R/integratedAnalysis5.R")
source( "~/ABARES/estimationCode/R/predictMethod1.R")
source( "~/ABARES/estimationCode/R/DoubleCountDataPrep1.R")
source( "~/ABARES/estimationCode/R/Pub_drawPosteriorSamps7.R")
source( "~/ABARES/estimationCode/R/helperFuns4.R")
source( "~/ABARES/estimationCode/R/isdmUtils4.R")

meshy1 <- makeMesh( ras=simDat$covarBrick[[1]], max.n=c(2000,350), dep.range=1)#, expans.mult=2.5)

fm.isdm1 <- isdm( POdat=as.data.frame( simDat$PO), AAdat=as.data.frame( simDat$AA), PA=as.data.frame( simDat$PA),
             covarBrick=simDat$covarBrick, mesh=meshy1,
             responseNames=c( PO="anything", AA="AA", PA="PA"), 
             sampleAreaNames=c( AA="transectArea", PA="transectArea"),
             distributionFormula=~0+Altitude+Temperature,
             biasFormulas=list( PO=~1+dist2City, AA=~1,  PA=~1),
             control=list( int.prec=0.01, other.prec=1,
                           calcICs=FALSE,
                           prior.range=c(0.05,0.1), prior.space.sigma=c( 2,0.0001),
                           coord.names=c("x","y"),
                           n.threads=8))
summary( fm.isdm1$mod)

fm.isdm1$preds <- predict.isdm( fit=fm.isdm1, covarRaster=simDat$covarBrick, S=500, DCaverage=TRUE)
popEsti <- PopEstimate( preds=fm.isdm1$preds, probs=c(0.025,0.975))




