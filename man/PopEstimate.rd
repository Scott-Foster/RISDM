\name{PopEstimate}
\alias{PopEstimate}
\title{Generates an estimate of population size from an isdm}
\description{This function estimates the number of individuals within the spatial domain which has been predicted. Please be aware that it may only be sensible under very specific scenarios: no preferential sampling (or adjustment for) and a sensible detectability estimate. The latter is only available in RISDM by using double count data.}
\usage{
 PopEstimate(preds, intercept.terms=NULL, control=NULL)
}
\arguments{
\item{preds}{An list as generated from \code{predict.isdm}.}
\item{intercept.terms}{The terms to use as part of the prediction process. If NULL (default) it is assumed that the intercepts have already been included as part of the prediciton process (see \code{predict.isdm}).}
\item{control}{Control arguments for how the population estimates are to be calculated. See details.}
}
\details{
 The function simply sums up, over the prediction raster, the cell-intensities. It does so for each posterior sample, thus giving a posterior summary.
 
 Two sets of predictions are given. The first three elements of the return object \emph{do not} incorporate the sampling/space distribution, whilst the second three have this extra source of variation included.

 Users should be aware that obtaining population estimates from these models needs a certain number of "stars to align". Just two of them are a sensible sampling plan, where the intercept terms are unbiassedly estimated, and for detectability a DC dataset should be used.
 
 The control parameters consist of
 * \code{probs} specifying the percentiles of the predictive distribution to return. Default is c(0.025,0.975) for 95\% confidence/prediction intervals.
 * \code{winsor} should each cell's MC samples be Winsorised prior to calculating the population estimate. Default is TRUE.
 * \code{tail} if Winsorising, what tail should be winsorised? Default is 'upper' and other options are 'lower' and 'both'.
 * \code{percet} if Winsorising, what proportion of data in the tail should be adjusted? Default is 0.01, corresponding to the cut-off values being 0.01 and 0.99 for lower and upper respectively.
}
\value{
 NULL
}

\seealso{\code{\link{isdm}}, \code{\link{predict.isdm}}}

\author{Scott D. Foster}

