
\name{plot.isdm}
\alias{plot.isdm}
\title{Plotting residuals for isdm objects}
\usage{
  \method{plot}{isdm}(x, covars, nFigRow=1, ask=TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"isdm"}, usually, a result of a call to \code{\link{isdm}}.}
  \item{covars}{a raster, or something that inheritits from a raster (raster*), to calculate residuals for PO data.}
  \item{nFigRow}{an integer indicating the number of rows for the figure. If less than the number of data types, there will be multiple figures. If greater, there will be empty rows. Default is 1, indicating that each data type gets its own page.}
  \item{ask}{a boolean, TRUE (default) if the user should be asked before plotting the residuals for the next data type. FALSE for no asking. Function always returns the state to not prompting on function exit.}
  \item{\dots}{argments passed to other sub-functions.  Currently not implemented.}
}
\description{
  This function is a method for class \code{isdm} objects.
}
\details{
  This method first calls the \code{\link{residuals.isdm}} method and then plots them in a (hopefully) useful manner. If the type of plot that you require is not part of the default, then please calculate the residuals (using the residuals.isdm method) and create the plot manually.
}
\seealso{
 \code{\link{isdm}}, \code{\link{predict.isdm}}, \code{\link{plot.isdm}}.
}

\author{Scott D. Foster}
\references{
  Dunn, P. K. and Smyth, G. K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics, \emph{5}, 236--244
}



