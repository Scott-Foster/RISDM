\name{PopEstimate}
\alias{PopEstimate}
\title{Generates an estimate of population size from an isdm}
\description{This function estimates the number of individuals within the spatial domain which has been predicted. Please be aware that it may only be sensible under very specific scenarios: no preferential sampling (or adjustment for) and a sensible detectability estimate. The latter is only available in RISDM by using double count data.}
\usage{
 PopEstimate(preds, probs=c(0.025,0.975), intercept.terms=NULL)
}
\arguments{
\item{preds}{An list as generated from \code{predict.isdm}.}
\item{probs}{The limits of the interval estimate returned. Expressed as quantiles.}
\item{intercept.terms}{The terms to use as part of the prediction process. If NULL (default) it is assumed that the intercepts have already been included as part of the prediciton process (see \code{predict.isdm}).}
}
\details{
 The function simply sums up, over the prediction raster, the cell-intensities. It does so for each posterior sample, thus giving a posterior summary.

 Users should be aware that obtaining population estimates from these models needs a certain number of "stars to align". Just two of them are a sensible sampling plan, where the intercept terms are unbiassedly estimated, and for detectability a DC dataset should be used.
}
\value{
 NULL
}

\seealso{\code{\link{isdm}}, \code{\link{predict.isdm}}}

\author{Scott D. Foster}

