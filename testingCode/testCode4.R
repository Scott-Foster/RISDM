
############################################
####  Checking the current development version of RISDM

#Complex boundary taken from victoria.  Project and extent used for simple too.

library( RISDM)

library( raster)
library( "RandomFields")
RFoptions(install="no")
source( "ReadData2.R")

rm( VicKoala_PO, VicKoala_PO_df, VicKoala_DC_df, VicKoala_SC_df, VicKoala_DC, VicKoala_SC, do.plots, tmp)

bdry.ras <- covars[[1]]
bdry.ras <- addLayer( bdry.ras, bdry.ras)
bdry.ras[[1]] <- boundaries( bdry.ras[[1]])
bdry.ras[[2]] <- boundaries( bdry.ras[[2]])
values( bdry.ras[[1]]) <- 1
values( bdry.ras[[1]])[!is.na( values( bdry.ras[[1]]))] <- 1
values( bdry.ras[[2]])[!is.na( values( bdry.ras[[2]]))] <- 1
names( bdry.ras) <- c("square","vic")
plot( bdry.ras)

bdry.ras <- aggregate( bdry.ras, 10, mean, na.rm=TRUE)
plot( bdry.ras)

#source( "~/ABARES/estimationCode/R/simulate2.R")

my.coefs <- list(dist=c(-1,1,-1), bias=c(-2,-1.5))

simDat1 <- simulate.isdm( expected.n.PO=3000, n.PA=200, n.AA=300, n.DC=200, 
                          coefs=my.coefs, transect.size=0.5,   #really huge transects -- not real just ismulation.
                          control=list( set.random.seed=TRUE, random.seed=787, do.plot=FALSE,
                                        sd=1, range=50000, addRandom=TRUE),
                          rasterBoundary=bdry.ras$vic)
crs( simDat1$covarBrick) <- crs( bdry.ras$square)

plot( simDat1$covarBrick)

dev.off()

meshy1 <- makeMesh( ras=simDat1$covarBrick[[1]], max.n=c(1250,250), dep.range=10000, expans.mult=4, offset=50000, expandRegion=TRUE)

#containers for the square and vic simulations.
fm1 <- list()

#with the PO data only using plugin estimates for now.

fm1[[1]] <- isdm( observationList=list( POdat=as.data.frame( simDat1$PO)),
                  covarBrick=simDat1$covarBrick, 
                  mesh=meshy1,
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
                                DCmethod="plugin"))#TalyorsLinApprox"))

fm1[[2]] <- isdm( observationList=list( DCdat=as.data.frame( simDat1$DC)),
                  covarBrick=simDat1$covarBrick, 
                  mesh=meshy1,
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
                                DCmethod="plugin"))#TalyorsLinApprox"))

fm1[[3]] <- isdm( observationList=list( AAdat=as.data.frame( simDat1$AA)),
                  covarBrick=simDat1$covarBrick, 
                  mesh=meshy1,
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
                                DCmethod="plugin"))#TalyorsLinApprox"))

fm1[[4]] <- isdm( observationList=list( PAdat=as.data.frame( simDat1$PA)),
                  covarBrick=simDat1$covarBrick, 
                  mesh=meshy1,
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
                                DCmethod="plugin"))#TalyorsLinApprox"))


fm1[[5]] <- isdm( observationList=list( POdat=as.data.frame( simDat1$PO), 
                                        DCdat=as.data.frame( simDat1$DC),
                                        AAdat=as.data.frame( simDat1$AA)),
#                                        PAdat=as.data.frame( simDat1$PA)),
                  covarBrick=simDat1$covarBrick, 
                  mesh=meshy1,
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
                                DCmethod="plugin"))#TalyorsLinApprox"))

summary( fm1[[1]]$mod)
summary( fm1[[2]]$mod)
summary( fm1[[3]]$mod)
summary( fm1[[4]]$mod)
summary( fm1[[5]]$mod)

tmp <- list()
for( ii in 1:length( fm1)){
  fm1[[ii]]$preds <- predict.isdm( fit=fm1[[ii]], covarRaster=simDat1$covarBrick, S=500, DCaverage=TRUE)
  tmp[[ii]] <- stack( crop( simDat1$covarBrick$Intensity, fm1[[ii]]$preds$mean.field), fm1[[ii]]$preds$mean.field)
  raster::plot( tmp[[ii]], nc=2)
  fm1[[ii]]$popEsti <- PopEstimate( preds=fm1[[ii]]$preds, probs=c(0.025,0.975))
  print( fm1[[ii]]$popEsti)
}



####checking PO stuff
tmp <- rasterize( x=simDat1$PO, y=bdry.ras$vic, background=0, fun='count')
tmp <- addLayer( tmp, raster( terra::cellSize( terra::rast( tmp))))
names( tmp) <- c( "count", "cellArea")

tmp1 <- SpatialPointsDataFrame( coords=coordinates( tmp), data=cbind( as.data.frame( tmp), as.data.frame( simDat1$covarBrick)))
simDat1$PO_as_AA <- as.data.frame( tmp1)

tmpTmp <- glm( count~1+Altitude+Temperature+dist2City, data=simDat1$PO_as_AA, family=poisson(link='log'), offset=log( cellArea))

summary( tmpTmp)
tmpPreds <- predict.glm( tmpTmp, newdata = simDat1$PO_as_AA, type='response')  #offset is added in when using the offset argument!?
tmpRaster <- rasterFromXYZ( cbind( simDat1$PO_as_AA[,c("x","y")], tmpPreds))
plot( tmpRaster)
sum( tmpPreds, na.rm=TRUE)

