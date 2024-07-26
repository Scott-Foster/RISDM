\name{simulateData.isdm}
\alias{simulateData.isdm}
\title{Simulates data of different types from an underlying point process.}
\description{This function simulates some data based on a number of fairly strict assumptions. }
\usage{
  simulateData.isdm( pop.size=10000,
			distForm=~-1+var1, biasForm=~1,
			Intercept=NULL, distCoefs=NULL, biasCoefs=NULL, DC.pi=NULL,
			n.PO=300, n.PA=150, n.AA=50, n.DC=50,
			rasterBoundary=NULL, covarBrick=NULL,
			transect.size=0.125,
			Intensity=NULL,
			control=list())
}

\arguments{
\item{pop.size}{The number of individuals expected throughout the sampling area. Ignored if Intercept is given (arguments play the same role).}
\item{distForm}{The formula used to describe the distribution model. It is important to not include an intercept here.}
\item{biasForm}{The formula used to describe the bias model.}
\item{Intercept}{The intercept to use in the distribution model. Default is NULL, which indicates that (by default) intercept is specified by the number of individuals.}
\item{distCoefs}{The coefficients to use in the model for the distribution. Note that this does not include an intercept.}
\item{biasCoefs}{The coefficients used to define the model for the bias.}
\item{DC.pi}{The detection probability for double count data. The same detection probability is assumed for each observer.}
\item{n.PA}{The number of presence-absence (PA) data.}
\item{n.AA}{The number of abundance (AA) data.}
\item{n.DC}{The number of double counted (DC) transects.}
\item{rasterBoundary}{SpatRaster, whose NA pattern defines the boundary of the simulation area. Ignored if rasterCovars is supplied.}
\item{covarBrick}{SpatRaster containing data described in distForm and biasForm.}
\item{transect.size}{The area covered by each of the PA, AA, and DC transects.}
\item{Intensity}{A user-supplied intensity surface to use instead of a model. If NULL (default) then the argument is ignored. If a raster (only other option) then it is used and distForm and distCoefs are ignored. Note however, that bias is still added to the PO sampling.}
\item{control}{A list of control arguments for the simulation. May be partially specified. See Details.}
}
\details{
 This function generates some fake data. It was written largely for internal testing purposes, but is made available just in case others find it useful (or useful to hack). The simulation of the random fields is facilitate by Skip Woolley's fftGPsim() function, which is copied (with permission) into the RISM package and then slightly refactored by Wen-Hsi Yang -- thanks Skip and Wen-Hsi!.
 
  The algorithm proceeds by first defining a grid, then calculating distances, Matern covariances, then decomposing covariances and simulating using decomposition. By far and away, the slowest piece is the decomposition -- a Cholesky is used.  Dense grids, greater than about 50x50, could be annoyingly slow.
 
 The control argument may contain the following elements:
\describe{
  \item{set.random.seed}{Should the random seed be set to random.seed? Default is FALSE: it should not be set.}
  \item{random.seed}{The value of the random seed. Ignored if set.random.seed==FALSE. Default is 787, a big plane.}
  \item{raster.dim}{The size of the raster to simulate over. Default is rep( 100,2), which corresponds to 100 cells in either direction.}
  \item{doPlot}{Should the elements arising from the simulation be plotted.  Default is TRUE, they should.}
  \item{addRandom}{Should a spatially explicit random effect be added to the intensity surface. Default is TRUE, it should.}
  \item{sd}{The standard deviation of the random effect.}
  \item{range}{The effective range of the random effect.}
}
}
\value{
 An object of class `simISDMdata', which is just a list with elements for each of the different data types, covariates, coefficients and so on. Names should be pretty self-explaining.
}

\seealso{\code{\link{isdm}}, \code{\link{predict.isdm}}}

\author{Scott D. Foster}

