\name{isdm}
\alias{isdm}
\title{Fits a spatial SDM to disparate data sources}
\description{ This function fits a SDM to presence-only and/or presence-absence and/or abundance-absence data amd/or double-count transect data. It does so using INLA (Rue et al.; 2009), a Poisson process approximateion, and is wrapped up in the integrated methods described in Fltetcher et al. (2019) and Isaac et al. (2020). Currently, data in any of the previously mentioned three data types can be used. The model rests on the assumption that the underlying distribution of the species is a log-Guassian Cox process (a Poisson process with a spatial random effect). The code started as an altered version of that supplied with Isaac et al. (2020), but has undergone method development and extension since then.}
\usage{
 isdm( observationList=list( POdat=NULL, PAdat=NULL, AAdat=NULL, DCdat=NULL),
	covars=NULL,
	habitatArea=NULL,
        mesh=NULL,
        responseNames=NULL,
        sampleAreaNames=NULL, 
        distributionFormula=NULL,
        biasFormula=NULL,
        artefactFormulas=list( PO=NULL, PA=~1, AA=~1, DC=~1),
        DCobserverInfo=list( SurveyID="surveyID", Obs1="Obs1", Obs2="Obs2", Both="Both"), #info for observer stuff
        control=list())
}
\arguments{
\item{ observationList}{A named list of data.frames containing the data from each type of data. For list element POdat (for presence-only data) this is a data.frame of occurance locations. For PAdat and AAdat (presence-absence and abundance data) these are data.frames containing the spatial location and any covariates associated with those particular observations (see artefactFormulas argument). For DCdat (double count data) the number observed from each observer and from both observers must be present (names in DCobserverInfo argument), and also relevant covariates.}
\item{ covars}{A RasterBrick (or RasterStack) object containing all the raster information (mostly covariates but see habitatArea argument) in the distribution and bias models.}
\item{ habitatArea}{A character giving the name of a layer of covars that corresponds to the amount of area within each cell of the raster brick that is suitable habitat for the species. This is an area value, units should match those that define the sampling area for the PA, AA, and DC data.}
\item{ mesh}{A mesh of class inla.mesh. Easiest to create using the makeMesh() function, which is a simplified wrapper to the inla.mesh.2d() function from the INLA package. Default is NULL, in which case this function will terminate with an error. Careful consideration is generally required for creating meshes. See ?makeMesh for some further details.}
\item{ responseNames}{A named character vector. The elements define the outcome variables in each of the datasets in POdat, PAdat, and AAdat. As an example, responseNames=c(PO="myPresences", AA="myAbundances") will pick out the variable myPresences from POdat and myAbundances from AAdat. The example assumes that there is no PA data. Note that no response name is required/requested for DCdata. For DCdata, the information is passed via the DCobserverInfo argument. For responseNames argument, the default is NULL, a value that will throw a (hopefully) helpful error.}
\item{ sampleAreaNames}{A named character vector giving the names of the area searched for each data type. This will be used as an offset value in each of the component models. No value is needed for PO data (they are points). Default is NULL, which will throw a (hopefully) helpful error.}
\item{ distributionFormula}{A formula that describes how the species relates to the covariates \emph{but without any potential bias}. Typically, this will include environmental covariates, basis-expanded versions of them (e.g. quadratics) and interactions. Terms in this formula must be present in covars. We note that basis expansions are useful and effective (e.g. poly() orthogonal polynomials), but these terms are standardised by default (away from being orthonormal).}
\item{biasFormula}{A formula describing the bias of the PO data. The default is an intercept-only formula, which encapsulates the assumption that the PO contain no observer-bias -- the pattern varies only with the distribution of the organisms. Terms in this formula must be present in covars. Like the distributionFormula, basis expansions are allowed, but are once again standardised.}
\item{artefactFormulas}{A list of formulas describing variation that only occurs within its respective data type. As a concrete example, AA data may be gathered from multiple institutes/sources and may all have different detection probabilities. In this case an entry of AA=~1+dataSource may suffice. All terms included in these formulas must be present in their respective data set. While models that use basis expansion will run, its precise meaning is unknown -- for now it is *not recommended* that basis expansions are used for any of the artefact formulas.}
\item{DCobserverInfo}{ A named list containing column names of DCdat. These columns contain: 1) a survey indicator (a different detection probability will be assumed for each observer in each survey); 2) the number of animals that observer 1 ("Obs1") saw, but observer 2 did not; 3) the number of animals that observer 2 ("Obs2") saw, but observer 1 did not, and; 3) the number that both observers saw ("Both").}
\item{control}{ A named list of arguments that control the estimation process. See Details below.}
}
\details{ This function is firmly based on the methods described in Isaac et al. (2020). It even started its learly ife as being based on their code-base. The code has grown substantially in terms of scope and methods, but the idea and intuition remains the same. One of the differences is that prediction is now performed using the INLA function \code{INLA::inla.posterior.samples}, as a separate analysis step. See \code{\link{predict.isdm}}.

All formulas may contain anything that a regular formula can parse (e.g. \code{\link{poly}} and \code{\link{bs}} and \code{\link{I}}) as well as interactions etc. However, please note that such terms may have unpredictable effects in the artefact formulae -- be careful (and probably do so by consturcting design matrices outside of this function). This warning should definitely be heeded when trying to include the same basis expansion for a common variable in both the distribution/bias and artefact formulae.

Control arguments for mdoel estimation come from the argument "control", which is a list. The elements of the control list are

\describe{
\item{n.threads}{The number of cores/threads to spread the computing over. Default is the number of cores available minus 1}
\item{coord.names}{The names for the spatial coordinates. Default is c("Easting","Northing")}
\item{prior.list}{A list of priors for each effect in the model. Has to be explicitly in the form required by \code{\link{control.fixed}} from the INLA package. No checking of this control parameter is performed. In particular, the coefficient labels have to be correct -- users may want to run a generic model first to obtain these. If prior.list is specified, then prior.mean, int.prec, other.prec are all ignored.}
\item{prior.mean}{The expectation for Gaussian prior for the `fixed' effects. Default value is zero. Ignored if prior.list is specified.}
\item{int.prec}{The precision for the Gaussian prior for the intercept value. Default value is 0.001. Ignored if prior.list is specified.}
\item{other.prec}{The precision for the Gaussian prior for the `fixed' effects. Default is 0.1.  Does not include any of the intercepts. Ignored if prior.list is specified.}
\item{calcICs}{Boolean value indicating if the performance measures should be calculated.  In particular, the marginal log-likelihood (model evidence).}
\item{prior.range}{The values to define the prior for the range parameter of the spatial random variable. Uses penalised complexity priors (Simpson et al. (2017). The default might be useful for a unit square and is c(0.2,0.1), which implies that the prior specifies that there is pr(practical range < 0.2) = 0.1. \emph{This will need to be changed for nearly every application.}}
\item{prior.space.sigma}{The values to define the prior for the variance of the spatial random variable. Uses penalised complexity priors (Simpson et al. (2017). The default is c( 1,0.01) meaning that pr( sigma>1)=0.01.  That is there isn't a whole lot of spatial variance. \emph{This may need to be changed for nearly every application.}}
\item{verbose}{Should INLA be run in verbose mode for debugging? Default is FALSE (no extraneous printing).}
\item{addRandom}{Should the random effect be included in the model. Default is TRUE (random effect included in the model).}
\item{standardiseCovariates}{Should the covariates for the *distribution* and *bias* formulae (not artefact formulae) be standardised before fitting. This can help with numerical issues, and it helps put covariate effects on a similar/same scale for specifying priors. Note that covariates for artefact formulae are not standardised and that the user should perform their own standardisation there, if needed.}
\item{returnStack}{Should the INLA data stack be returned as part of the model output. Default is TRUE -- it is returned (at the expense of memory). Note that this object is needed for use in \code{\link{plot.isdm}}.}
\item{DCmethod}{Which method should be used to estimate the detectability of each observer. Default is "TaylorsLinApprox", which approximates the model. If this doesn't work well, then the other option is "plugin", which estimates the probabilities prior to estimation.}
}
}
\value{An object of class "isdm". This is just a list with elements that contain: 
\item{mod}{The INLA model.}
\item{distributionFormula}{The formula for the distribution model.}
\item{biasFormula}{The formula for the bias in the PO data.}
\item{artefactFormulas}{The list of formulas for the artefact terms for each of the data types.}
\item{mesh}{The INLA mesh object that is used to compute the random spatial field.}
\item{control}{The isdm control parameters used for the model's fit}
\item{responseNames}{The names of the response data for the different data types. Mostly used internally.}
\item{data}{A list of two elements. The first is a raster brick with the variables used in the model fit (expanded as per model.matrix). The second element is itself a list and contains the data for the artefact models (again expanded as per model.matrix).}
}

\seealso{ \code{\link{makeMesh}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}, \code{\link{predict.isdm}}, \code{\link{plot.isdm}}}

\author{Scott D. Foster}
\references{
  Rue, H., Martino, S. and Chopin, N. (2009) Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society: Series B (Statistical Methodology) \emph{71}, 319--392.

  Fletcher Jr., R. J., Hefley, T. J., Robertson, E. P., Zuckerberg, B., McCleery, R. A. and Dorazio, R. M. (2019) A practical guide for combining data to model species distributions Ecology, \emph{100}, e02710

  Isaac, N. J., Jarzyna, M. A., Keil, P., Dambly, L. I., Boersch-Supan, P. H., Browning, E., Freeman, S. N., Golding, N., Guillera-Arroita, G., Henrys, P. A., Jarvis, S., Lahoz-Monfort, J., Pagel, J., Pescott, O. L., Schmucki, R., Simmonds, E. G. and O’Hara, R. B. (2020) Data Integration for Large-Scale Models of Species Distributions Trends in Ecology & Evolution, \emph{35}, 56--67

  Simpson, D. P., Rue, H., Riebler, A., Martins, T. G. and Sørbye, S. H. (2017) Penalising model component complexity: A principled, practical approach to constructing priors Statistical Science, \emph{32}, 1--28
}
