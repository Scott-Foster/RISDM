
\name{plot.isdm}
\alias{plot.isdm}
\title{Plotting residuals for isdm objects}
\usage{
  \method{plot}{isdm}(x, covarRaster, \dots)
}
\arguments{
  \item{x}{an object of class \code{"isdm"}, usually, a result of a call to \code{\link{isdm}}.}
  \item{covarRaster}{a raster, or something that inheritits from a raster (raster*), to calculate residuals for PO data.}
  \item{ask}{a boolean, TRUE (default) if the user should be asked before plotting the residuals for the next data type. FALSE for no asking. Function always returns the state to not prompting on function exit.}
  \item{\dots}{argments passed to other sub-functions.  Currently only implemented is S from \code{\link{predict.isdm}}, specifying the number of draws to sample from the posterior distribution.}
}
\description{
  This function is a \code{\link{methods}} for class \code{isdm} objects.
}
\details{
  Residuals for DC, AA and PA data are directly applications of the randomised quantile residual (RQR) approach of Dunn and Smyth (1996). The residuals are taken from the posterior expectation of the fitted model (taken from INLA directly). RQR residuals are general, but in our experience (annecdotal evidence) may not be all-that-sensitive to departures from model assumptions.  That is, the model may not be as good as the residual plots suggest.

  The residuals for the PO data are formed spatially. The area is gridded into a raster, whose number of PO observations defines the value for the cell. Predictions are then made (using \code{\link{predict.isdm}}) and the RQR calculated.

  For all types of residuals, no attempt has been made to account for (e.g.) spatial autocorrelation, amongst other known departures. These residuals are intended to be a guide rather than a definitive solution to diagnostics.

  Note that the \code{\link{predict.isdm}} method is used, rather than taking the prediction from the INLA object, as the predict method will correctly scale to the cell area.  
}
\seealso{
 \code{\link{isdm}}, \code{\link{predict.isdm}}.
}

\author{Scott D. Foster}
\references{
  Dunn, P. K. and Smyth, G. K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics, \emph{5}, 236--244
}



