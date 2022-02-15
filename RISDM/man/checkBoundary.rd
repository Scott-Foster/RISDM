\name{checkMesh}
\alias{checkMesh}
\title{Performs elementary checks for a mesh.}
\description{ This is a helper function. It is designed to simplify the process of creating a good mesh to perform spatial computing over. Hence, it is useful prior to calling \code{\link{isdm}} and in conjunction with \code{\link{makeMesh}}. The function takes a mesh and provides some primitive diagnostics to see if its component triangles are likely to be suitable. We hope that it simplifies the mesh creation process. These checks are performed only on the inner domain -- the outer domain will often (?) have quite different sized triangles.}
\usage{
 checkMesh( mesh, hull)
}
\arguments{
\item{mesh}{A mesh generated from \code{\link{makeMesh}} or \code{\link{inla.mesh.2d}}.}
\item{hull}{A spatial boundary corresponding to the boundary of the mesh's inner and outer domains. The easiest way to obtain this object is through using \code{\link{makeMesh}} or \code{\link{inla.nonconvex.hull}}. If using makeMesh, then the necessary object is stored in the list element 'hull'}
}
\details{
 This function does three things: 1) produces a plot of the mesh (duplicating what you probably already have done), 2) produces a histogram of the distribution of the areas of the mesh's component triangles, and 3) produces a histogram of the distribution of all the internal angles for all the triangles. The rule of thumbs, mostly from Krainsky et al. (2017), suggest that there should be
 
 1) No areas where triangles are commonly larger -- check the plot of the mesh.
 2) No overly large triangles -- check the distribution of their areas.
 3) All triangles should be relatively equalateral, with no triangles having overly small or overly large internal angles -- check the distribution of the angles.
 
 While these diagnostics **may** be useful, we are aware that they might also **not be** useful. Also there may be other more useful diagnostics. Our advice, if time and computing allow, is to use a number of different meshes and see what sensitivity the results have.
}
\value{
 NULL
}

\seealso{\code{\link{makeBoundary}}, \code{\link{isdm}}, \code{\link{checkMesh}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}}

\author{Scott D. Foster}

\references{
  Bakka, Haakon 2017, Mesh Creation including Coastlines, accessed 20 January 2022, <https://haakonbakkagit.github.io/btopic104.html#3_The_simplest_mesh>.
  
  Krainski, E.; Gómez-Rubio, V.; Bakka, H.; Lenzi, A.; Castro-Camilo, D.; Simpson, D.; Lindgren, F. & Rue, H. Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA Chapman & Hall/CRC Press, 2019
  
  Lindgren, F.; Rue, H. & Lindström, J. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach Journal of the Royal Statistical Society: Series B (Statistical Methodology), Blackwell Publishing Ltd, 2011, 73, 423-498
}
