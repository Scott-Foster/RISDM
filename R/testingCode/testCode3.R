
source( "~/ABARES/estimationCode/R/simulate2.R")

my.coefs <- list(dist=c(-1,-1,1), bias=c(-2,-0.5))

##########################################
####  Is the simulation process OK??

#Abundance (AA) checking using a glm
simDat <- simulate.isdm( n.AA=500, coefs=my.coefs, control=list( set.random.seed=TRUE, random.seed=767, do.plot=TRUE,
                                                       sd=1.5, range=4, addRandom=TRUE), rasterBoundary=bdry.ras$vic)
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
                                                                 sd=1.5, range=4, addRandom=TRUE), rasterBoundary=bdry.ras$vic)

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
source( "~/ABARES/estimationCode/R/simulate2.R")

my.coefs <- list(dist=c(-1,1,-1), bias=c(-2,-1.5))

simDat1 <- simulate.isdm( expected.n.PO=1000, n.PA=250, n.AA=250, n.DC=250, 
                          coefs=my.coefs, transect.size=0.25, 
                          control=list( set.random.seed=TRUE, random.seed=121, do.plot=TRUE,
                                        sd=1, range=50000, addRandom=TRUE),
                          rasterBoundary=bdry.ras$square)
simDat2 <- simulate.isdm( expected.n.PO=1000, n.PA=250, n.AA=250, n.DC=250, 
                         coefs=my.coefs, transect.size=0.25, 
                         control=list( set.random.seed=TRUE, random.seed=121, do.plot=TRUE,
                                       sd=1, range=50000, addRandom=TRUE),
                         rasterBoundary=bdry.ras$vic)

plot( simDat1$covarBrick)
plot( simDat2$covarBrick)

dev.off()
##using new code
source( "~/ABARES/estimationCode/R/integratedAnalysis7.R")
source( "~/ABARES/estimationCode/R/predictMethod1.R")
source( "~/ABARES/estimationCode/R/DoubleCountDataPrep2.R")
source( "~/ABARES/estimationCode/R/Pub_drawPosteriorSamps7.R")
source( "~/ABARES/estimationCode/R/helperFuns4.R")
source( "~/ABARES/estimationCode/R/isdmUtils6.R")
source( "~/ABARES/estimationCode/R/bookMeshDual.R")

meshy1 <- makeMesh( ras=simDat1$covarBrick[[1]], max.n=c(1500,400), dep.range=10000, expans.mult=4, offset=50000)
meshy2 <- makeMesh( ras=simDat2$covarBrick[[1]], max.n=c(1500,400), dep.range=10000, expans.mult=4, offset=50000)

#containers for the square and vic simulations.
fm1 <- fm2 <- list()

#with the DC data only

#doesn't want to fit -- too few DC observations.
fm1[[1]] <- isdm( observationList=list( DCdat=as.data.frame( simDat1$DC)),
                  covarBrick=simDat1$covarBrick, 
                  mesh=meshy1,
                  responseNames=c( DC="blah"),
                  sampleAreaNames=c( DC="transectArea"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=NULL,
                  artefactFormulas=list( DC=~1+Survey),
                  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE, 
                                DCmethod="TalyorsLinApprox"))

fm2[[1]] <- isdm( observationList=list( DCdat=as.data.frame( simDat2$DC)),
                  covarBrick=simDat2$covarBrick, 
                  mesh=meshy2,
                  responseNames=c( DC="blah"),
                  sampleAreaNames=c( DC="transectArea"),
                  distributionFormula=~0+Altitude+Temperature,
                  biasFormula=NULL,
                  artefactFormulas=list( DC=~1+Survey),
                  DCobserverInfo=list( SurveyID="Survey", Obs1="Obs1", Obs2="Obs2", Both="Both"),
                  control=list( int.prec=0.01, other.prec=1,
                                calcICs=FALSE,
                                prior.range=c(1000,0.1), prior.space.sigma=c( 2,0.1),
                                coord.names=c("x","y"),
                                n.threads=8,
                                addRandom=TRUE, 
                                DCmethod="TalyorsLinApprox"))

summary( fm1[[1]]$mod)
summary( fm2[[1]]$mod)


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

## with the PO data only

fm.isdm4 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO)), 
  #AAdat=as.data.frame( simDat$AA),
#  PAdat=as.data.frame( simDat$PA)), 
  #DCdat=as.data.frame( simDat$DC),
  covarBrick=simDat$covarBrick, mesh=meshy1,
  responseNames=c( PO="PO"),#PO="anything", AA="AA", PA="PA"),
  sampleAreaNames=c( PO="transectArea"),#AA="transectArea", PA="transectArea", DC="transectArea"),
  distributionFormula=~0+Altitude+Temperature+dist2City,
  biasFormula=~-1,
  artefactFormulas=list( PO=~-1),#PO=~1+dist2City, AA=~1,  PA=~1, DC=~1),
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
print( fm.isdm4$mod$summary.hyperpar)

print( fm.isdm1a$mod$summary.fixed)
print( fm.isdm1b$mod$summary.fixed)
print( fm.isdm2$mod$summary.fixed)
print( fm.isdm3$mod$summary.fixed)
print( fm.isdm4$mod$summary.fixed)

####  That all looks pretty good for this one simulation.  ALthough the hyperparams for the PO are a little out compared to the rest. 
####  Lots of variation in places though.

###  I've ahd previous simulations where the 'fancier' DC approach doesn't work.
###  The fixed effects for many of the models cover the true value, although sometimes only just...

# PO + AA
fm.isdm5 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
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
fm.isdm6 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
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
fm.isdm7 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
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

fm.isdm8 <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
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

fm.isdm8a <- isdm( observationList=list( POdat=as.data.frame( simDat$PO), 
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
                                DCmethod="TaylorsLinApprox"))

print( fm.isdm2$mod$summary.hyperpar) #AA only model
print( fm.isdm5$mod$summary.hyperpar)
print( fm.isdm6$mod$summary.hyperpar)
print( fm.isdm7$mod$summary.hyperpar)
print( fm.isdm8$mod$summary.hyperpar)
print( fm.isdm8a$mod$summary.hyperpar)

print( fm.isdm2$mod$summary.fixed)  #AA only model
print( fm.isdm5$mod$summary.fixed)
print( fm.isdm6$mod$summary.fixed)
print( fm.isdm7$mod$summary.fixed)
print( fm.isdm8$mod$summary.fixed)
print( fm.isdm8a$mod$summary.fixed)

### prediction....
fm.isdm8$preds <- predict.isdm( fit=fm.isdm8, covarRaster=simDat$covarBrick, S=500, DCaverage=TRUE)
raster::plot( fm.isdm8$preds$mean.field)
popEsti8 <- PopEstimate( preds=fm.isdm8$preds, probs=c(0.025,0.975))
print( popEsti8)
print( popEsti8$mean / simDat$expected.pop.size)

