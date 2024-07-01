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
library( terra)

## ----readCovars, eval=TRUE,fig.cap="Environmental covariate data for the gamba grass example. Accessibility (ACC) measuring effective distance for humans to travel. Digital elevation model (DEM) giving the altitude. Soil moisture at the root zone (SMRZ) measuring the wetness where it matters for plants. See text for details and references."----
filenames <- paste0( system.file("extdata", package="RISDM"), 
                   c("/GambaExample_sqrtACC_23Mar28.tif", 
                     "/GambaExample_sqrtDEM_23Mar28.tif", 
                     "/GambaExample_SMRZ_23Mar28.tif"))
covars <- rast( filenames)
names( covars) <- c("ACC","DEM","SMRZ")
print( colnames( crds( covars)))  #to see how things are labelled.
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
checkMesh( my.mesh)

## ----makeBadMesh1, eval=TRUE, fig.height=3------------------------------------
my.mesh.bad <- makeMesh( covars$DEM, max.n=c(250, 30), dep.range=3, doPlot=FALSE)
checkMesh( my.mesh.bad)

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

