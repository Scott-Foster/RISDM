\name{isdm}
\alias{isdm}
\title{Fits a spatial SDM to disparite data sources}
\description{ This function fits a SDM to presence-only and/or presence-absence and/or abundance-absence data amd/or double-count transect data. It does so using INLA (Rue et al.; 2009), a Poisson process approximateion, and is wrapped up in the integrated methods described in Fltetcher et al. (2019) and Isaac et al. (2020). Currently, data in any of the previously mentioned three data types can be used. The model rests on the assumption that the underlying distribution of the species is a log-Guassian Cox process (a Poisson process with a spatial random effect). The code started as an altered version of that supplied with Isaac et al. (2020), but has undergone method development and extension since then.}
\usage{
 isdm( observationList=list( POdat=NULL, PAdat=NULL, AAdat=NULL, DCdat=NULL),
	covarBrick=NULL,
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
\item{ covarBrick}{A RasterBrick (or RasterStack) object containing all the covariates in the distribution model.}
\item{ mesh}{A mesh of class inla.mesh. Easiest to create using the makeMesh() function, which is a simplified wrapper to the inla.mesh.2d() function from the INLA package. Default is NULL, in which case this function will terminate with an error. Careful consideration is generally required for creating meshes. See ?makeMesh for some further details.}
\item{ responseNames}{A named character vector. The elements define the outcome variables in each of the datasets in POdat, PAdat, and AAdat. As an example, responseNames=c(PO="myPresences", AA="myAbundances") will pick out the variable myPresences from POdat and myAbundances from AAdat. The example assumes that there is no PA data. Note that no response name is required/requested for DCdata. For DCdata, the information is passed via the DCobserverInfo argument. For responseNames argument, the default is NULL, a value that will throw a (hopefully) helpful error.}
\item{ sampleAreaNames}{A named character vector giving the names of the area searched for each data type. This will be used as an offset value in each of the component models. No value is needed for PO data (they are points). Default is NULL, which will throw a (hopefully) helpful error.}
\item{ distributionFormula}{A formula that describes how the species relates to the covariates \emph{but without any potential bias}. Typically, this will include environmental covariates, basis-expanded versions of them (e.g. quadratics) and interactions. Terms in this formula must be present in covarBrick.}
\item{biasFormula}{A formula describing the bias of the PO data. The default is an intercept-only formula, which encapsulates the assumption that the PO contain no observer-bias -- the pattern varies only with the distribution of the organisms. Terms in this formula must be present in covarBrick.}
\item{artefactFormulas}{A list of formulas describing variation that only occurs within its respective data type. As a concrete example, AA data may be gathered from multiple institutes/sources and may all have different detection probabilities. In this case an entry of AA=~1+dataSource may suffice. All terms included in these formulas must be present in their respective data set.}
\item{DCobserverInfo}{ A named list containing column names of DCdat. These columns contain: 1) a survey indicator (a different detection probability will be assumed for each observer in each survey); 2) the number of animals that observer 1 ("Obs1") saw, but observer 2 did not; 3) the number of animals that observer 2 ("Obs2") saw, but observer 1 did not, and; 3) the number that both observers saw ("Both").}
\item{control}{ A named list of arguments that control the estimation process. See Details below.}
}
\details{ This function is firmly based on the methods described in Isaac et al. (2020). It even started its life as being based on their code-base. The code has grown in terms of scope and methods, but the idea and intuition remains the same. 

One of the differences is that prediction is now performed using the INLA function \code{INLA::inla.posterior.samples}, as a separate analysis step. See \code{\link{predict.isdm}}.

Control arguments for mdoel estimation come from the argument "control", which is a listl. The elements of the control list are

\describe{
\item{n.threads}{The number of cores/threads to spread the computing over. Default is the number of cores available minus 1}
\item{coord.names}{The names for the spatial coordinates. Default is c("Easting","Northing")}
\item{prior.mean}{The expectation for Gaussian prior for the `fixed' effects. Default value is zero.}
\item{int.prec}{The precision for the Gaussian prior for the intercept value. Default value is 0.0001.}
\item{other.prec}{The precision for the Gaussian prior for the `fixed' effects. Default is 0.0001.  Does not include the intercept.}
\item{calcICs}{Boolean value indicating if the performance measures should be calculated.  In particular, the marginal log-likelihood (model evidence).}
\item{prior.range}{The values to define the prior for the range parameter of the spatial random variable. Uses penalised complexity priors (Simpson et al. (2017). The default might be useful for a unit square and is c(0.2,0.1), which implies that the prior specifies that there is pr(practical range < 0.2) = 0.1. \emph{This will need to be changed for nearly every application.}}
\item{prior.space.sigma}{The values to define the prior for the variance of the spatial random variable. Uses penalised complexity priors (Simpson et al. (2017). The default is c( 1,0.01) meaning that pr( sigma>1)=0.01.  That is there isn't a whole lot of spatial variance. \emph{This may need to be changed for nearly every application.}}
\item{verbose}{Should INLA be run in verbose mode for debugging? Default is FALSE (no extraneous printing).}
\item{addRandom}{Should the random effect be included in the model. Default is TRUE (random effect included in the model).}
\item{returnStack}{Should the INLA data stack be returned as part of the model output. Default is FALSE -- it is not returned to save memory.}
\item{DCmethod}{Which method should be used to estimate the detectability of each observer. Default is "TaylorsLinApprox", which approximates the model. If this doesn't work well, then the other option is "plugin", which estimates the probabilities prior to estimation.}
}
}
\value{An object of class "isdm". This is just a list with elements that contain: 
\item{mod}{The INLA model.}
\item{distributionFormula}{The formula for the distribution model.}
\item{biasFormula}{The formula for the bias in the PO data.}
\item{artefactFormulas}{The list of formulas for the artefact terms for each of the data types.}
\item{mesh}{The INLA mesh object that is used to compute the random spatial field.}
}

\seealso{ \code{\link{makeMesh}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}, \code{\link{predict.isdm}}}

\author{Scott D. Foster}
\references{
  Rue, H., Martino, S. and Chopin, N. (2009) Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society: Series B (Statistical Methodology) \emph{71}, 319--392.

  Fletcher Jr., R. J., Hefley, T. J., Robertson, E. P., Zuckerberg, B., McCleery, R. A. and Dorazio, R. M. (2019) A practical guide for combining data to model species distributions Ecology, \emph{100}, e02710

  Isaac, N. J., Jarzyna, M. A., Keil, P., Dambly, L. I., Boersch-Supan, P. H., Browning, E., Freeman, S. N., Golding, N., Guillera-Arroita, G., Henrys, P. A., Jarvis, S., Lahoz-Monfort, J., Pagel, J., Pescott, O. L., Schmucki, R., Simmonds, E. G. and O’Hara, R. B. (2020) Data Integration for Large-Scale Models of Species Distributions Trends in Ecology & Evolution, \emph{35}, 56--67

  Simpson, D. P., Rue, H., Riebler, A., Martins, T. G. and Sørbye, S. H. (2017) Penalising model component complexity: A principled, practical approach to constructing priors Statistical Science, \emph{32}, 1--28
}
