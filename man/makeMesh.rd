\name{makeMesh}
\alias{makeMesh}
\title{Creates a mesh to fit an \code{\link{isdm}} over.}
\description{ This is a helper function. It is designed to simplify the creation of a mesh to perform spatial computing over. Hence, it is useful prior to calling \code{\link{isdm}}. The function takes a raster and, using the INLA package in conjunction to some judicious defaults, finds a mesh governed by the arguments. We hope that it simplifies mesh creation, in that it attempts to find `an acceptable' mesh.}
\usage{
 makeMesh( ras, max.n=NULL, dep.range=NULL, expandRegion=TRUE, expans.mult=NULL, hull.res=100, max.edge=NULL, cutoff=NULL, offset=NULL, doPlot=TRUE, ...)
}
\arguments{
\item{ras}{A raster object whose patterns of non-NA values define the analysis domain and whose NA values define the non-domain. The actual contents of the raster do not matter and will be ignored (except whether cells contain an NA or not).}
\item{max.n}{An integer vector of 2 elements. It describes the maximum number of mesh nodes to place in the inner domain and the outer domain. The default (c500,200) is chosen quite arbitrarily -- users will want to consider this carefully.}
\item{dep.range}{Technically: The assumed effective range of the spatial random effect in the SDM. Here the effective range is the distance required to reduce correlation to 0.1. Intuitively: The distance required before sampling produces a more-or-less independent observation. Note that this is range is not inclusive of the covariate effects. Default (when argument is NULL) is taken to be 1/3 of the diagonal length of the raster's extent. This value is based on the spatial properties of the raster rather than any property of the species data themselves. As such, the default may be poor. It may be that a preliminary analysis, using \code{\link{isdm}}, could be used to obtain a reasonable value for this argument.}
\item{expandRegion}{DE-ACTIVATED. If you don't want expansion then set expans.mult to be some small positive number. Should the survey region be expanded beyond the boundary defined by the raster? Default is TRUE, and should nearly always be TRUE. Honestly, this parameter and the functionality is a legacy of trying bug-hunt.}
\item{expans.mult}{The inner domain is defined by a nonconvex hull expansion around the non-NA values in ras. The expansion amount can be linked to the dep.range and a useful default could be 1.5 (Bakka). This gives the expansion to be 1.5*dep.range.}
\item{hull.res}{The resolution of the nonconvex hull. This is largely a technical argument to make sure that computation works. Default is 100, which you should only adjust if a warning message suggests that it is inadequate. Smaller expans.mults are likely to need larger hull.res.}
\item{max.edge}{A numeric vector of 2 elements. The elements describe the largest allowable edge within the inner domain and the outer domain. The default is NULL, in which case the values of c(0.2,0.5)*dep.range are used. Note that for good results, the max edges in the inner domain should not be large with respect to the dep.range (Bakka).}
\item{cutoff}{Scalar numeric. This scalar describes the mesh building algorithms starting smallest distance between any two nodes (inner and outer). If NULL (default) a value of 0.2*max.edge is used.}
\item{offset}{The amount to expand the nonconvex hull defining the inner domain. This amount defines the outer boundary of the outer domain. Default is to extend by the dep.range. Note that this argument is an absolute value, not a multiple (even though the default is dep.range).}
\item{ doPlot}{Boolean indicating whether the mesh should be plotted. We highly recommend that all meshes are, at some point, to make sure that the mesh is behaving as you think it should be.}
\item{ ...}{Additional arguments for, and passed directly to, \code{\link{inla.mesh.2d}}. There is no checking of these arguments, but they shouldn't include those listed above.}
}
\details{
 This function relies heavily on the INLA package. In particular, nearly all of the heavy lifting is done by the \code{\link{inla.mesh.2d}} function. Users should look at \code{?inla.mesh.2d} for more details about that process.

 The creation of the outer area is important to ensure that the random spatial effect within the model used by \code{\link{isdm}} is properly specified and is less likely to exhibit edge effects.

 The number of nodes should be chosen carefully: too many nodes and the computational burden within isdm will be too great. Too few and the approximation used by isdm, through \code{\link{inla}} will deviate from reality.

 Specifying an acceptable mesh is an artform, somewhat akin to black magic it seems. To aleviate some of the resulting pain, we use rules-of-thumb from Bakka (2017) and Krianski et al. (2019). These rules motivate our defaults, but they should always be corefully considered and the resulting mesh be carefully critiqued. To help with the critique see \cite{\link{checkMesh}}.
 
 Note that the arguments are non-independent; sometimes tweaking one will do nothing until another argument it changed. This is a direct result of the underlying function \code{\link{inla.mesh.2d}} being quite smart and sensible and uses these values as limits, not goals. This should be entirely expected given the arguments' names.
}
\value{
 An extended inla.mesh object describing the mesh on which to perform computing within the region of interest. The extension has two extra list elements: 1) the boundary to the non-NA entries in the raster; and 2) the nonconvex hull around the non-NA elements.
}

\seealso{\code{\link{isdm}}, \code{\link{checkMesh}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}}

\author{Scott D. Foster}

\references{
  Bakka, Haakon 2017, Mesh Creation including Coastlines, accessed 20 January 2022, <https://haakonbakkagit.github.io/btopic104.html#3_The_simplest_mesh>.
  
  Krainski, E.; Gómez-Rubio, V.; Bakka, H.; Lenzi, A.; Castro-Camilo, D.; Simpson, D.; Lindgren, F. & Rue, H. Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA Chapman & Hall/CRC Press, 2019
  
  Lindgren, F.; Rue, H. & Lindström, J. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach Journal of the Royal Statistical Society: Series B (Statistical Methodology), Blackwell Publishing Ltd, 2011, 73, 423-498
}
