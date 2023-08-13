\name{simulateData.isdm}
\alias{simulateData.isdm}
\title{Simulates data of different types from an underlying point process.}
\description{This function simulates some data based on a number of fairly strict assumptions. }
\usage{
 simulateData.isdm( expected.pop.size=400000, expected.n.PO=300, n.PA=150, n.AA=50, n.DC=50,
                           coefs=list(dist=c(-1.5,-0.25,0.75), bias=c(-2,-1.5)), 
                           DC.pis=matrix( c(0.8,0.76, 0.7,0.73, 0.82,0.67), nrow=3, ncol=2, byrow=TRUE),
                           transect.size = 0.125, #a proportion of cell size.
                           rasterBoundary=NULL,
                           control=list())
}
\arguments{
\item{expected.pop.size}{The number of individuals expected throughout the sampling area.}
\item{expected.n.PO}{The number of observed Presence-only (PO) observations in the sampling area.}
\item{n.PA}{The number of presence-absence (PA) data.}
\item{n.AA}{The number of abundance (AA) data.}
\item{n.DC}{The number of double counted (DC) transects.}
\item{coefs}{The coefficients for the distribution model: intercept and each of the two covariates.}
\item{bias}{The coefficients for the bias model: intercept and the bias covariate.}
\item{DC.pis}{The detection probability for each of the two DC observers.}
\item{transect.size}{The area covered by each of the PA, AA, and DC transects.}
\item{rasterBoundary}{Raster, whose NA pattern defines the boundary of the simulation area.}
\item{control}{A list of control arguments for the simulation. May be partially specified. See Details.}
}
\details{
 This function generates some fake data. It was written largely for internal testing purposes, but is made available just in case others find it useful (or useful to hack). It is pretty rudimentary. The simulation of the random fields is facilitate by Skip Woolley's fftGPsim() function, which is copied (with permission) into the RISM package -- thanks Skip.
 
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
 An object of class `simISDMdata', which is just a list with elements for each of the different data types, covariates and the expected population size.
}

\seealso{\code{\link{isdm}}, \code{\link{predict.isdm}}}

\author{Scott D. Foster}

