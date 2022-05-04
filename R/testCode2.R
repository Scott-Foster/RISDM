
source( "~/ABARES/estimationCode/R/simulate1.R")

my.coefs <- list(dist=c(-1,-1,2), bias=c(-2,-2.5))

##########################################
####  Is the simulation process OK??

#Abundance (AA) checking using a glm
simDat <- simulate.isdm( n.AA=500, coefs=my.coefs, control=list( set.random.seed=TRUE, random.seed=767, do.plot=FALSE,
                                                       sd=1.5, range=4, addRandom=TRUE))
predDat <- raster::as.data.frame( simDat$covarBrick)
predDat$transectArea <- prod( raster::res( simDat$covarBrick))

tmpDat <- raster::as.data.frame( simDat$AA)
tmpDat <- cbind( tmpDat, raster::extract( simDat$covarBrick, tmpDat[,c("x","y")]))
fm.AA <- glm( AA~1+Altitude+Temperature, family=poisson(), offset=log(transectArea), data=tmpDat)
summary( fm.AA)  #seems to be pretty good!  Upto spatial RE confounding?
pred.AA <- predict( fm.AA, newdata=predDat, type='response')
print( sum( pred.AA) / simDat$expected.pop.size)

#Presence-absence (PA) checking using a glm
simDat <- simulate.isdm( n.PA=500, coefs=my.coefs, control=list( set.random.seed=TRUE, random.seed=767, do.plot=FALSE,
                                                                 sd=1.5, range=4, addRandom=TRUE))

tmpDat <- as.data.frame( simDat$PA)
tmpDat <- cbind( tmpDat, raster::extract( simDat$covarBrick, tmpDat[,c("x","y")]))
fm.PA <- glm( PA~1+Altitude+Temperature, family=binomial( link="cloglog"), offset=log( transectArea), data=tmpDat)
summary( fm.PA)  #seems to be pretty good!  Upto spatial RE confounding?
pred.PA <- exp( predict( fm.PA, newdata=predDat, type='link'))
print( sum( pred.AA) / simDat$expected.pop.size)

####  Simulation seems to be appropriate.  Although this doesn't check the PO data (as it is biassed and a simple model like above won't work)

############################################
####  Checking the current development version of RISDM

#library( RISDM)  #should get overwritten by the sourcing of following files
library( raster)
source( "~/ABARES/estimationCode/R/simulate1.R")

my.coefs <- list(dist=c(-1,-0.5,1), bias=c(-2,-1))

simDat <- simulate.isdm( expected.n.PO=1000, n.PA=250, n.AA=250, n.DC=250, 
                         coefs=my.coefs, transect.size=0.01, 
                         control=list( set.random.seed=TRUE, random.seed=12, do.plot=TRUE,
                                       sd=0.25, range=1.5, addRandom=TRUE))
plot( simDat$covarBrick)

dev.off()
##using new code
source( "~/ABARES/estimationCode/R/integratedAnalysis7.R")
source( "~/ABARES/estimationCode/R/predictMethod1.R")
source( "~/ABARES/estimationCode/R/DoubleCountDataPrep2.R")
source( "~/ABARES/estimationCode/R/Pub_drawPosteriorSamps7.R")
source( "~/ABARES/estimationCode/R/helperFuns4.R")
source( "~/ABARES/estimationCode/R/isdmUtils5.R")

meshy1 <- makeMesh( ras=simDat$covarBrick[[1]], max.n=c(500,200), dep.range=0.5, expans.mult=2, offset=1)

#with the DC data only

fm.isdm1a <- isdm( observationList=list( #POdat=as.data.frame( simDat$PO), 
                  #AAdat=as.data.frame( simDat$AA),
                  #PA=as.data.frame( simDat$PA), 
                  DCdat=as.data.frame( simDat$DC)),
                  covarBrick=simDat$covarBrick, 
                  mesh=meshy1,
                  responseNames=c( DC="blah"),#PO="anything", AA="AA", PA="PA"),
                  sampleAreaNames=c( DC="transectArea"),#AA="transectArea", PA="transectArea", DC="transectArea"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=NULL,
                  artefactFormulas=list( DC=~1+Survey),#PO=~1+dist2City, AA=~1,  PA=~1, DC=~1),
                  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE, 
                                DCmethod="TalyorsLinApprox"))

#with the DC data and plug in estimates for observers (internal plug-in values).

fm.isdm1b <- isdm( observationList=list(  #POdat=as.data.frame( simDat$PO), 
  DCdat=as.data.frame( simDat$DC)),
  #PA=as.data.frame( simDat$PA), 
  #DCdat=as.data.frame( simDat$DC),
  covarBrick=simDat$covarBrick, 
  mesh=meshy1,
  responseNames=c( DC="blah"),#PO="anything", AA="AA", PA="PA"),
  sampleAreaNames=c( DC="transectArea"),#AA="transectArea", PA="transectArea", DC="transectArea"),
  distributionFormula=~0+Altitude+Temperature,
  biasFormula=NULL,
  artefactFormulas=list( DC=~-1+Survey),#PO=~1+dist2City, AA=~1,  PA=~1, DC=~1),
  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
  control=list( int.prec=0.01, other.prec=1,
                calcICs=FALSE,
                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                coord.names=c("x","y"),
                n.threads=8,
                addRandom=TRUE, 
                DCmethod="plugin"))

#with the AA data only

fm.isdm2 <- isdm( observationList=list( #POdat=as.data.frame( simDat$PO), 
                                        AAdat=as.data.frame( simDat$AA)),
                                        #PA=as.data.frame( simDat$PA), 
                                        #DCdat=as.data.frame( simDat$DC),
  covarBrick=simDat$covarBrick, 
  mesh=meshy1,
  responseNames=c( AA="AA"),#PO="anything", AA="AA", PA="PA"),
  sampleAreaNames=c( AA="transectArea"),#AA="transectArea", PA="transectArea", DC="transectArea"),
  distributionFormula=~0+Altitude+Temperature,
  biasFormula=NULL,
  artefactFormulas=list( AA=~1),#PO=~1+dist2City, AA=~1,  PA=~1, DC=~1),
#  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
  control=list( int.prec=0.01, other.prec=1,
                calcICs=FALSE,
                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                coord.names=c("x","y"),
                n.threads=8,
                addRandom=TRUE, 
                DCmethod="plugin"))

## with the PA data only

fm.isdm3 <- isdm( observationList=list( #POdat=as.data.frame( simDat$PO), 
                                        #AAdat=as.data.frame( simDat$AA),
                                        PAdat=as.data.frame( simDat$PA)), 
                                        #DCdat=as.data.frame( simDat$DC),
  covarBrick=simDat$covarBrick, mesh=meshy1,
  responseNames=c( PA="PA"),#PO="anything", AA="AA", PA="PA"),
  sampleAreaNames=c( PA="transectArea"),#AA="transectArea", PA="transectArea", DC="transectArea"),
  distributionFormula=~0+Altitude+Temperature,
  biasFormula=NULL,
  artefactFormulas=list( PA=~1),#PO=~1+dist2City, AA=~1,  PA=~1, DC=~1),
  #  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
  control=list( int.prec=0.01, other.prec=1,
                calcICs=FALSE,
                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                coord.names=c("x","y"),
                n.threads=8,
                addRandom=TRUE))

#sd=1.5, range=4
print( fm.isdm1a$mod$summary.hyperpar)
print( fm.isdm1b$mod$summary.hyperpar)
print( fm.isdm2$mod$summary.hyperpar)
print( fm.isdm3$mod$summary.hyperpar)

print( fm.isdm1a$mod$summary.fixed)
print( fm.isdm1b$mod$summary.fixed)
print( fm.isdm2$mod$summary.fixed)
print( fm.isdm3$mod$summary.fixed)

###  So for these data the fancier DC approach doesn't work (1a versus 1b hyperparams when there is spatial noise).
###  I don't know why.  Drop fancy approach for now.
###  Curiously, the fixed effects for many of the models cover the true value, although sometimes only just...

# PO + AA
fm.isdm4 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
                                        AAdat=as.data.frame( simDat$AA)),
                                        #PA=as.data.frame( simDat$PA), 
                                        #DCdat=as.data.frame( simDat$DC),
  covarBrick=simDat$covarBrick, 
  mesh=meshy1,
  responseNames=c( AA="AA"),# PA="PA"),
  sampleAreaNames=c( AA="transectArea"),# PA="transectArea", DC="transectArea"),
  distributionFormula=~0+Altitude+Temperature,
  biasFormula=~1+dist2City,
  artefactFormulas=list( AA=~1),#,  PA=~1, DC=~1),
  #  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
  control=list( int.prec=0.01, other.prec=1,
                calcICs=FALSE,
                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                coord.names=c("x","y"),
                n.threads=8,
                addRandom=TRUE))

# PO + DC
fm.isdm5 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
                                        DCdat=as.data.frame( simDat$DC)),
                  #PA=as.data.frame( simDat$PA), 
                  #DCdat=as.data.frame( simDat$DC),
                  covarBrick=simDat$covarBrick, 
                  mesh=meshy1,
                  responseNames=c(DC="blah"),#c( AA="AA"),# PA="PA"),
                  sampleAreaNames=c( DC="transectArea"),# PA="transectArea", DC="transectArea"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=~1+dist2City,
                  artefactFormulas=list( DC=~1+Survey),#,  PA=~1, DC=~1),
                  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE))
# PO + PA
fm.isdm6 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
                                        PAdat=as.data.frame( simDat$PA)),
                  #PA=as.data.frame( simDat$PA), 
                  #DCdat=as.data.frame( simDat$DC),
                  covarBrick=simDat$covarBrick, 
                  mesh=meshy1,
                  responseNames=c( PA="PA"),# PA="PA"),
                  sampleAreaNames=c( PA="transectArea"),# PA="transectArea", DC="transectArea"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=~1+dist2City,
                  artefactFormulas=list( PA=~1),#,  PA=~1, DC=~1),
                  #  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE))


## PO + AA + DC + PA too.

fm.isdm7 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
                  AAdat=as.data.frame( simDat$AA),
                  PAdat=as.data.frame( simDat$PA), 
                  DCdat=as.data.frame( simDat$DC)),
                  covarBrick=simDat$covarBrick, 
                  mesh=meshy1,
                  responseNames=c( AA="AA", PA="PA"),
                  sampleAreaNames=c( AA="transectArea", PA="transectArea", DC="transectArea"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=~1+dist2City,
                  artefactFormula=list( AA=~1,  PA=~1, DC=~1+Survey),
                  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),#sorting out observer stuff, offset and coef
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(0.1,0.01), prior.space.sigma=c( 10,0.001),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE, 
                                DCmethod="plugin"))

print( fm.isdm2$mod$summary.hyperpar)
print( fm.isdm4$mod$summary.hyperpar)
print( fm.isdm5$mod$summary.hyperpar)
print( fm.isdm6$mod$summary.hyperpar)
print( fm.isdm7$mod$summary.hyperpar)

print( fm.isdm2$mod$summary.fixed)
print( fm.isdm4$mod$summary.fixed)
print( fm.isdm5$mod$summary.fixed)
print( fm.isdm6$mod$summary.fixed)
print( fm.isdm7$mod$summary.fixed)





### prediction....
fm.isdm7$preds <- predict.isdm( fit=fm.isdm7, covarRaster=simDat$covarBrick, S=500, DCaverage=TRUE)
raster::plot( fm.isdm7$preds$mean.field)
popEsti7 <- PopEstimate( preds=fm.isdm7$preds, probs=c(0.025,0.975))
print( popEsti7)
print( popEsti7$mean / simDat$expected.pop.size)

