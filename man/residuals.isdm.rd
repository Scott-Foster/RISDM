
\name{residuals.isdm}
\alias{residuals.isdm}
\title{Calculating Randomised Quantile Residuals for isdm objects}
\usage{
  \method{residuals}{isdm}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{"isdm"}, usually, a result of a call to \code{\link{isdm}}.}
  \item{\dots}{argments passed to other sub-functions.  Currently not implemented.}
}
\description{
  This function is a method for class \code{isdm} objects.
}
\details{
  Residuals for DC, AA and PA data are directly applications of the randomised quantile residual (RQR) approach of Dunn and Smyth (1996). The residuals are taken from the posterior expectation of the fitted model (taken from INLA directly). RQR residuals are general, but in our experience (annecdotal evidence) may not be all-that-sensitive to departures from model assumptions.  That is, the model may not be as good as the residual plots suggest.

  The residuals for the PO data are formed spatially. The area is gridded into a raster, whose number of PO observations defines the value for the cell. Predictions are then made (using \code{\link{predict.isdm}}) and the RQR calculated.

  For all types of residuals, no attempt has been made to account for (e.g.) spatial autocorrelation, amongst other known departures. These residuals are intended to be a guide rather than a definitive solution to diagnostics.

  Note that the \code{\link{predict.isdm}} method is used, rather than taking the prediction from the INLA object, as the predict method will correctly scale to the cell area.  
}
\value{An object of class "isdm". This is just a list with elements that contain: 
  \item{DC}{A data.frame containing the point predictions, the observation and the RQR residual.}
  \item{AA}{As per DC but for AA data.}
  \item{PA}{As per DC but for PA data.}
  \item{PA}{A list of two elements. The first (ras) is a raster containing the RQR residuals organised into spatial locations. The second (POresids) is a data.frame with fitted values, observations and RQR.}
}

\seealso{
 \code{\link{isdm}}, \code{\link{predict.isdm}}, \code{\link{plot.isdm}}.
}

\author{Scott D. Foster}
\references{
  Dunn, P. K. and Smyth, G. K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics, \emph{5}, 236--244
}



