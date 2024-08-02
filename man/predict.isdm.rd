\name{predict.isdm}
\alias{predict.isdm}
\title{Makes predictions from an isdm object}
\description{ This function predicts from an isdm object from the \code{\link{isdm}} function.  It does so by sampling from the approximate posterior of the model and produces a posterior raster.}
\usage{
 \method{predict}{isdm}( object, covars, habitatArea=NULL, S=500, intercept.terms=NULL, 
			n.threads=NULL, n.batches=1, includeRandom=TRUE, includeFixed=TRUE, includeBias=FALSE, 
			type="intensity", confidence.level=0.95, quick=FALSE, ...)
}

\arguments{
\item{object}{An object of class isdm, as obtained from \code{isdm}. No default.}
\item{covars}{A rasterBrick (or rasterStack) object containing all the covariates in the distribution model.}
\item{habitatArea}{A character giving the name of a layer of \code{covars} that corresponds to the amount of area within each cell of the raster brick that is suitable habitat for the species. Typically, useage here will correspond with use in the \code{\link{isdm}} function too.}
\item{S}{The number of posterior samples to take. Default is 500 samples, which is likely to be small for serious applications.}
\item{intercept.terms}{Vector of strings indicating which terms in the model should be included as intercepts. An example might be c("AA_Intercept","AA_Intercept:surveyIDdonna") meaning that the coefficient for AA_Intercept and for the interaction AA_Intercept:surveyIDdonna will both be added to each of the predictions. If NULL (default), the function will not add any intercept. To ensure that you get the correct intercept name, it is easiest to take the text for the correct term from fit$mod$names.fixed or from the \code{\link{summary}} method.}
\item{n.threads}{How many threads to spread the computation over. Default is NULL, where the number used to estimate the model (arugment "fit") is used.}
\item{n.batches}{How many batches should the prediction be split into. This is largely a hack to overcome memory issues when needing to use a large S. The number of posterior samples will be broken into n.batches groups. Avoids an error, I'd suggest not using this unless you see the error.}
\item{includeRandom}{Should the random spatial effect be included in the predictions? Default is TRUE, as it nearly always should be (unless you are trying to understand the contribution of the terms).}
\item{includeFixed}{Should the fixed effects, and/or which ones, be invluded in the predictions? Default is TRUE indicating that the fixed effects (and intercepts from the intercept.terms argument) are included in the predictions. If FALSE, then no fixed effects nor intercepts are included. An alternative option is to specify the name of (a) variable(s) in the formula (and also a name of a layer in the covars rasterStack). Such specification allows for partial effects to be predicted.}
\item{includeBias}{Should the sampling bias be included in the predictions? Default is FALSE, it is not included. This term is nearly always not-interesting in terms of figuring out what is where. However, it could be interesting to see where the search effort has been placed. Please be aware that includeBias=TRUE will force the intercept.PO to be added to the linear predictor. So, please do not include a non-NULL intercept.terms argument as that will make multiple intercepts to be included in the prediction.}
\item{type}{The type (scale) of prediction. Choices are "intensity" for the parameter of the log-Guass Cox process, "probability" for the probability of having any 1 observation in the prediction cell, or "link" for the linear predictor.}
\item{confidence.level}{The intended coverage for the confidence intervals. Default is 0.95 for 95 percent intervals (ranging from 0.025 percentile to 0.975 percentile).}
\item{quick}{Should the predictions be performed using MC sampling of the posterior (quick=FALSE and the default) or should a quick and dirty value be obtained using the mean of each parameters' posterior. See details for more.}
\item{...}{Not implemented}
}

\details{ This function is a isdm specific interface to \code{INLA::inla.posterior.samples}. The function generates samples, selects which ones should be included for predicting and then performs the necessary machinations to do the predictions. All predictions are for the grid of covarRaster, including the area of each cell. That is the prediction is for the number of individuals (from the point process) within a cell.

The covariate data, within \code{covars}, is treated to the same preparation as the covariate data for the model fit using \code{\link{isdm}}. In particular, the internal standardisation in \code{\link{isdm}} is used in \code{predict.isdm} too -- by subtracting the grid's mean and dividing by it's standard deviation. This means that the use of the *same* rasters in both sets will have the desired results. However, if a different raster is provided for prediction than at estimation, then results may be unpredictable. This includes prediction into new geographical areas and/or prediction at a different resolution. If these kinds of predictions are required, then it is recommended that the user supplies their own raster layers for estimation and predictions which have been carefully constructed to avoid potential numerical under or overflow. Briefly, this means that the values of the cells shouldn't be overly large or small. These supplied raster layers should have the same scaling and the user should set control$standardiseCovariates=FALSE in the \code{\link{isdm}} estimation process. The standardisation in the prediction process will follow that used in the estimation process.}

Whilst the \code{quick} argument allows for some visualisation of the model, it is not complete nor is it likely to be hugely accurate. For accuracy that you'll need quick=FALSE. This feature is included only to enable model development for large models.

\value{A list containing the following elements:
\item{field}{The predictions and summaries thereof. The summaries are for the cell-wise posterior: median, lower interval, upper interval, mean, and standard deviation. Note that, in general, the median is a more robust measure of central tendancy than the mean.}
\item{cell.samples}{All the S posterior draws for each cell.}
\item{fixedSamples}{All the S posterior draws of the fixed effects.}
\item{fixed.names}{The names of the fixed effects in the samples, fixedSamples}
\item{predLocats}{The locations where the predictions are made}
\item{confidence.limits}{The percentiles for the bounds of the confidence limits.}
}

\seealso{ \code{\link{makeMesh}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}, \code{INLA::inla.posterior.samples}}

\author{Scott D. Foster}
