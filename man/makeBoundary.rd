\name{makeBoundary}
\alias{makeBoundary}
\title{Creates a spatial polygon surrounding the non-NA region(s) of a raster layer.}
\description{ This is a helper function. It is designed to simplify the creation of a boundary SpatialPolygon that is suitable for use in mesh creation, and hence for \code{\link{isdm}}. The function takes a raster layer and, using the sp package, finds the edge locations. These locations are assembled into a polygon, whose internal `holes' have been removed. If told, it can also do some (desireable) thinning of the raster so that the resulting polygon is not too jaggered and information rich.  This thinning will speed up computing and make mesh generation easier.}
\usage{
 makeBoundary( r, lowerRes=TRUE, doPlot=TRUE)
}
\arguments{
\item{ r}{A raster layer containing NAs for all cells outside of the region of interest. No default, must be specified.}
\item{ lowerRes}{A boolean indicating whether the raster r should have its resolution lowered. Default is TRUE -- resolution will be down-graded. The amount of down-grading is done so that there are approximately 100 cells in each direction within the raster. Lowering the resolution can/will help producing computationally possible models based on interpretable meshes.}
\item{ doPlot}{Should the resulting polygon be plotted. Default (TRUE) indicates that it should.}
}
\details{
 This function relies heavily on the raster and sp packages. In particular, most of the heavy lifting is done by the \code{\link{rasterToPolygons}} and the \code{\link{aggregate}} functions.
}
\value{
 A SpatialPolygons object describing the boundary of the interesting region within the raster.
}

\seealso{\code{\link{makeMesh}}, \code{\link{isdm}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}}

\author{Scott D. Foster}
