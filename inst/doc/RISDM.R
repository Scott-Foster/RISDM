## ----prelim, echo = FALSE, results="hide"-------------------------------------
library( knitr)
opts_chunk$set(cache=FALSE, message = FALSE, comment = "", dev="pdf",
                      dpi=50, fig.show = "hold", fig.align = "center")

## ----setup1, eval=FALSE-------------------------------------------------------
#  install.packages("INLA",repos=c(getOption("repos"),
#                  INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#  #vignette not built to save computation
#  devtools::install_github( repo="Scott-Foster/RISDM", build_vignettes=FALSE)

## ----setup2-------------------------------------------------------------------
library( RISDM)
library( raster)

## ----readCovars, eval=TRUE,fig.cap="Environmental covariate data for the gamba grass example. Accessibility (ACC) measuring effective distance for humans to travel. Digital elevation model (DEM) giving the altitude. Soil moisture at the root zone (SMRZ) measuring the wetness where it matters for plants. See text for details and references."----
filenames <- paste0( system.file("extdata", package="RISDM"), 
                   c("/GambaExample_sqrtACC_23Mar28.tif", 
                     "/GambaExample_sqrtDEM_23Mar28.tif", 
                     "/GambaExample_SMRZ_23Mar28.tif"))
covars <- stack( filenames)
names( covars) <- c("ACC","DEM","SMRZ")
print( colnames( coordinates( covars)))  #to see how things are labelled.
plot( covars)

## ----readObs, eval=TRUE-------------------------------------------------------
gamba_PO <- readRDS( system.file("extdata", "Gamba_PO_23Mar28.RDS", 
                                   package="RISDM"))
gamba_PA <- readRDS( system.file("extdata", "Gamba_PA_23Mar28.RDS", 
                                 package="RISDM"))

str( gamba_PO)
str( gamba_PA)

## ----makeMesh1, eval=TRUE-----------------------------------------------------
my.mesh <- makeMesh( covars$DEM, max.n=c(1000, 350), dep.range=3, doPlot=FALSE)

## ----checkMesh1, eval=TRUE, fig.height=3--------------------------------------
checkMesh( my.mesh, my.mesh$hull, ras=covars$DEM)

## ----makeBadMesh1, eval=TRUE, fig.height=3------------------------------------
my.mesh.bad <- makeMesh( covars$DEM, max.n=c(250, 30), dep.range=3, doPlot=FALSE)
checkMesh( my.mesh.bad, my.mesh$hull, ras=covars$DEM)

## ----specifyDist, eval=TRUE---------------------------------------------------
#linear in SMRZ only
my.form <- ~0+SMRZ
#interacting between SMRZ and DEM
my.form <- ~0+SMRZ*DEM
#a B-spline regression basis with low degrees of freedom
my.form <- ~0+splines::bs(SMRZ, df=3)
#a(n orthogonal) polynomial in SMRZ and DEM
my.form <- ~0+poly(SMRZ,2)+poly(DEM,2)
#adding accessibility too (as per paper).
my.form <- ~0+poly(SMRZ,2)+poly(DEM,2)+acc

## ----biasForm, eval=TRUE------------------------------------------------------
#intercept only == no heterogeneity in search effort
my.biasForm <- 1
#regression spline in acc
my.biasForm <- ~1+splines::bs( acc, df=3)
#linear in acc
my.biasForm <- ~1+acc

## ----arteForm, eval=TRUE------------------------------------------------------
#observer differences as a sampling artefact
list( PA=~1+observer)
#intercept only == no sampling artefact 
#               == no heterogeneity within each data type
list( PA=~1)

## ----initFit, eval=TRUE-------------------------------------------------------
fm <- isdm( observationList=list( POdat=gamba_PO, 
                                  PAdat=gamba_PA),
            covars=covars, 
            mesh=my.mesh,
            responseNames=c( PO=NULL, PA="PA"),
            sampleAreaNames=c( PO=NULL, PA="Area"),
            distributionFormula=~0+poly( DEM, 2) + poly( SMRZ,2) + ACC,
            biasFormula=~1+ACC,
            artefactFormulas=list( PA=~1),
            control=list( coord.names=c("x","y")))

## ----contrEG, eval=TRUE-------------------------------------------------------
#the very vague everything
my.control <- list( prior.mean=0, int.sd=1000, other.sd=1000)
#(much) tighter prior for effects
my.control <- list( prior.mean=0, int.sd=1000, other.sd=0.1)
#the default
my.control <- list( prior.mean=0, int.sd=1000, other.sd=10)

## ----fos24Fit, eval=TRUE------------------------------------------------------
fm <- isdm( observationList=list( POdat=gamba_PO, 
                                  PAdat=gamba_PA),
            covars=covars, 
            mesh=my.mesh,
            responseNames=c( PO=NULL, PA="PA"),
            sampleAreaNames=c( PO=NULL, PA="Area"),
            distributionFormula=~0+poly( DEM, 2) + poly( SMRZ,2) + ACC,
            biasFormula=~1+ACC,
            artefactFormulas=list( PA=~1),
            control=list( coord.names=c("x","y"), 
                          int.sd=1000, other.sd=10, prior.mean=0,
                          prior.range=c(1,0.1), prior.space.sigma=c( 5,0.1)))

## ----noRandFit, eval=TRUE-----------------------------------------------------
fm.noRand <- isdm( observationList=list( POdat=gamba_PO, 
                                         PAdat=gamba_PA),
                   covars=covars, 
                   mesh=my.mesh,
                   responseNames=c( PO=NULL, PA="PA"),
                   sampleAreaNames=c( PO=NULL, PA="Area"),
                   distributionFormula=~0+poly( DEM, 2) + poly( SMRZ,2) + ACC,
                   biasFormula=~1+ACC,
                   artefactFormulas=list( PA=~1),
                   control=list( coord.names=c("x","y"), 
                                  int.sd=1000, other.sd=10, prior.mean=0,
                                  addRandom=FALSE))

## ----summ,eval=TRUE-----------------------------------------------------------
summary( fm)

## ----residPlots, eval=TRUE, fig.height=3.5, fig.cap="Residual plots for the gamba grass data, and the model containing quadratic effects and a spatial term."----
plot( fm, covars=covars, nFigRow=2, ask=FALSE)

## ----residPlots2, eval=TRUE, fig.height=3.5, fig.cap="Residual plots for the gamba grass data, and the model containing quadratic effects \textit{but no spatial term}."----
plot( fm.noRand, covars=covars, nFigRow=2, ask=FALSE)

## ----pred,eval=TRUE,fig.cap="Predictions from the quadratic model including random effects. Upper two rows are for the intensity, and bottom two rows are predictions of the probability of presence (within a raster cell)", fig.height=3.5----
#You should use a much(!) larger value of S.
#You may want to choose a larger value of n.threads too.
fm$preds <- predict( fm, covars=covars, S=50, 
                       intercept.terms="PA_Intercept")
plot( fm$preds$field)
  
#Predicting probability too
fm$preds.probs <- predict( fm, covars=covars, 
                           S=50, intercept.terms="PA_Intercept", 
                           type="probability")
plot( fm$preds.probs$field)

## ----interp, eval=TRUE, fig.height=3.5, fig.cap="Relationship with Soil Moisture (SMRZ). Black solid line is the median relationship and grey shaded area is the 95 percent CI."----
#the data for interpretation
#adding a temporary cell area layer
covarsForInter <- addLayer( covars, covars[[1]])
names( covarsForInter) <- c( names( covars), "tmp.habiArea")
#area is now constant with log(1)=0 contribution
values( covarsForInter$tmp.habiArea) <- 1 

#You could use a much(!) larger value of S.
interpPreds <- predict( fm, covars=covarsForInter, 
                        habitatArea="tmp.habiArea", S=50, 
                        includeFixed="SMRZ", includeRandom=FALSE, type="link")

#compile covariate and prediction
pred.df <- as.data.frame( cbind( SMRZ=values( covars$SMRZ), 
     values( interpPreds$field[[c("mu.median","mu.lower","mu.upper")]])))
#plot
pred.df <- pred.df[!is.na( pred.df$SMRZ),]
pred.df <- pred.df[order( pred.df$SMRZ),]
matplot( pred.df[,1], pred.df[,2:4], pch="", xlab="SMRZ", ylab="Effect", 
                                                main="Effect plot for SMRZ")
polygon( x=c( pred.df$SMRZ, rev( pred.df$SMRZ)), 
            c(pred.df$mu.upper, rev(pred.df$mu.lower)), 
                                                    col=grey(0.95), bor=NA)
lines( pred.df[,c("SMRZ","mu.median")], type='l', lwd=2)

## ----singleDataPO, eval=TRUE, fig.height=3.5, fig.cap="Intensity model and predictions from estimation using only PA data.", fig.height=3.5----
#PO data only
fm.PO <- isdm( observationList=list( POdat=gamba_PO),
            covars=covars, 
            mesh=my.mesh,
            responseNames=NULL,
            sampleAreaNames=NULL,
            distributionFormula=~0+poly( DEM, 2) + poly( SMRZ,2) + ACC,
            biasFormula=~1+ACC,
            artefactFormulas=NULL,
            control=list( coord.names=c("x","y"), 
                          int.sd=1000, other.sd=10, prior.mean=0,
                          prior.range=c(1,0.1), prior.space.sigma=c( 5,0.1)))
fm.PO$preds <- predict( fm.PO, covars=covars, S=50, 
                       intercept.terms="PO_Intercept")
plot( fm.PO$preds$field)

## ----singleDataPA, eval=TRUE, fig.height=3.5, fig.cap="Intensity model and predictions from estimation using only PA data.", fig.height=3.5----
#PA data only
fm.PA <- isdm( observationList=list( PAdat=gamba_PA),
            covars=covars, 
            mesh=my.mesh,
            responseNames=c( PA="PA"),
            sampleAreaNames=c( PA="Area"),
            distributionFormula=~0+poly( DEM, 2) + poly( SMRZ,2) + ACC,
            artefactFormulas=list( PA=~1),
            control=list( coord.names=c("x","y"), 
                          int.sd=1000, other.sd=10, prior.mean=0,
                          prior.range=c(1,0.1), prior.space.sigma=c( 5,0.1)))
fm.PA$preds <- predict( fm.PA, covars=covars, S=50, 
                       intercept.terms="PA_Intercept")
plot( fm.PA$preds$field)

## ----Tidy, eval=FALSE---------------------------------------------------------
#  #You may wish to tidy your workspace.
#  rm( covars, fm, fm.noRand, fm.PA, fm.PO, gamba_PA, gamba_PO, filenames,
#               my.biasForm, my.form, my.control, my.mesh, my.mesh.bad)

## ----sessionInfo, results = "asis", echo = FALSE------------------------------
toLatex(sessionInfo())

