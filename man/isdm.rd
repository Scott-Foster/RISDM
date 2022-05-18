\name{isdm}
\alias{isdm}
\title{Fits a spatial SDM to disparite data sources}
\description{ This function fits a SDM to presence-only and/or presence-absence and/or abundance-absence data. It does so using INLA (Rue et al.; 2009), the Poisson process method in Simpson et al. (2016), and wrapped up in the integrated methods described in Fltetcher et al. (2019) and Isaac et al. (2020). Currently, data in any of the previously mentioned three data types can be used. The model rests on the assumption that the underlying distribution of the species is a log-Guassian Cox process (a Poisson process with a spatial random effect). The code is an altered version of that supplied with Isaac et al. (2020).}
\usage{
 isdm( POdat=NULL, PAdat=NULL, AAdat=NULL, DCdat=NULL,
        covarBrick=NULL,
        mesh=NULL,
        boundary=NULL,
        responseNames=NULL,
        sampleAreaNames=NULL,
        distributionFormula=NULL,
        biasFormulas=list( PO=~1, PO=~1, AA=~1, DC=~1),  #an intercept for each data type.
        DCobserverInfo=list( SurveyID="surveyID", Obs1="Obs1", Obs2="Obs2", Both="Both"), #info for observer stuff
        control=list())
}
\arguments{
\item{ POdat}{A data.frame containing at least the spatial locations of the presences in the presence-only data. Default is NULL, implying that there is no PO data. Note that at least one of POdat, PAdat, AAdat and DCdat must be non-NULL.}
\item{ PAdat}{A data.frame containing at least the spatial locations and the presence-absence status at sampled locations. Default is NULL. See POdat.}
\item{ AAdat}{A data.frame containing at least the spatial locations and the abundance measure at sampled locations. Abundance measure must be able to include zero. Default is NULL. See POdat.}
\item{ DCdat}{A data.frame containing at least the spatial locations and the counts (plural) at sampled locations. Counts must have zeros measured (if observed). Needs three columns of counts: for the count observered by both observers, those only by observer 1, and those only by observer 2.}
\item{ covarBrick}{A RasterBrick (or RasterStack) object containing all the covariates used in the model.}
\item{ mesh}{A mesh of class inla.mesh. Easiest to create using the inla.mesh.2d() function from the INLA package, or a simplified wrapper for it from this package called makeMesh(). Default is NULL, in which case this function will terminate with an error. Careful consideration is generally required for creating meshes. See ?makeMesh for some further details.}
\item{ boundary}{A SpatialPolygon specifying the limit of the modelling area. PO, PA, and AA points should all be within this boundary. If NULL (default), the boudnary is guessed by taking the extent of the covarBrick object. It many situations this is unlikely to be what is wanted.}
\item{ responseNames}{A named character vector of length 1<=length(responseNames)<=3. The elements define the outcome variables in each of the datasets in POdat, PAdat, and AAdat. As an example, responseNames=c(PO="myPresences", AA="myAbundances") will pick out the variable myPresences from POdat and myAbundances from AAdat. The example assumes that there is no PA data (it would fall over otherwise). Note that no response name is required/requested for DCdata. For DCdata, the information is passed via the DCobserverInfo argument. For responseNames argument, the default is NULL, a value that will throw a (hopefully) helpful error.}
\item{ sampleAreaNames}{ A named character vector of length 1<=length(sampleAreaNames)<=3 giving the names of the area searched for each data type. This will be used as an offset value in each of the component models. No value is needed for PO data (they are points). Default is NULL, which will throw a (hopefully) helpful error.}
\item{ distributionFormula}{ A formula that describes how the species relates to the covariates \emph{but without any potential bias}. Typically, this will include environmental covariates, basis-expanded versions of them (e.g. quadratics) and interactions.}
\item{biasFormulas}{ A named list of formulas describing the sampling artefacts for each data type. The default is a list with intercept-only formulas, which encapsulates the assumption that the data types have a parallel version of the point process (on link scale). Including other terms can, for example, account for variation between different sampling methods.}
\item{DCoververInfo}{ A named list containing column names of DCdat. These columns contain: 1) a survey indicator (a different detection probability will be assumed for each observer in each survey); 2) the number of animals that observer 1 ("Obs1") saw, but observer 2 did not; 3) the number of animals that observer 2 ("Obs2") saw, but observer 1 did not, and; 3) the number that both observers saw ("Both").}
\item{control}{ A named list of arguments that control the estimation process. See Details below.}
}
\details{ This function is firmly based on the methods described in Isaac et al. (2020). It is even based on their code-base. There are a few changes and expansions to their original code, but the idea and intuition remains the same. Notably, double count data can now be utilised.

Control arguments for mdoel estimation come from the argument "control", which is a listl. The elements of the control list are

\describe{
\item{n.threads}{The number of cores/threads to spread the computing over. Default is the number of cores available minus 1}
\item{tag.pred}{The tag for the INLA model that defines the predictions. Default is "pred"}
\item{spat.index}{The name for the spatial term. Default is "i"}
\item{coord.names}{The names for the spatial coordinates. Default is c("Easting","Northing")}
\item{Meshpars}{A named list to be used to generate the mesh to perform computing over. Default is list( max.edge = 0.1, offset = 0, cutoff = 0.001, max.n=1000) but these should be adjusted depending on spatial scales and CRS. See ?INLA::inla.mesh.2d for a description of what the values do.}
\item{prior.mean}{The expectation for Gaussian prior for the `fixed' effects. Default value is zero.}
\item{int.prec}{The precision for the Gaussian prior for the intercept value. Default value is 0.0001.}
\item{other.prec}{The precision for the Gaussian prior for the `fixed' effects. Default is 0.0001.  Does not include the intercept.}
\item{calcICs}{Boolean value indicating if the performance measures should be calculated.  In particular, the marginal log-likelihood (model evidence).}
\item{prior.range}{The values to define the prior for the range parameter of the spatial random variable. Uses penalised complexity priors (Simpson et al. (2017). The default might be useful for a unit square and is c(0.2,0.1), which implies that the prior specifies that there is pr(practical range < 0.2) = 0.1. \emph{This will need to be changed for nearly every application.}}
\item{prior.space.sigma}{The values to define the prior for the variance of the spatial random variable. Uses penalised complexity priors (Simpson et al. (2017). The default is c( 1,0.01) meaning that pr( sigma>1)=0.01.  That is there isn't a whole lot of spatial variance. \emph{This may need to be changed for nearly every application.}}
\item{verbose}{Should INLA be run in verbose mode for debugging? Default is FALSE (no extraneous printing).}
}
}
\value{A list of two named elements, mod and preds. mod is an INLA model that contains all the information from the inla estimation. The second element, preds, is a RasterBrick containing the posterior expectations and standard deviations for the prediction locations.
}

\seealso{\code{\link{makeBoundary}}, \code{\link{makeMesh}}, \code{\link{inla.mesh.2d}}, \code{\link{inla}}}

\author{Scott D. Foster but based on code from Isaac et al. 2020}
\references{
  Rue, H., Martino, S. and Chopin, N. (2009) Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society: Series B (Statistical Methodology) \emph{71}, 319--392.

  Simpson, D., Illian, J. B., Lindgren, F., Sørbye, S. H. and Rue, H. (2016) Going off grid: computationally efficient inference for log-Gaussian Cox processes. Biometrika

  Fletcher Jr., R. J., Hefley, T. J., Robertson, E. P., Zuckerberg, B., McCleery, R. A. and Dorazio, R. M. (2019) A practical guide for combining data to model species distributions Ecology, \emph{100}, e02710

  Isaac, N. J., Jarzyna, M. A., Keil, P., Dambly, L. I., Boersch-Supan, P. H., Browning, E., Freeman, S. N., Golding, N., Guillera-Arroita, G., Henrys, P. A., Jarvis, S., Lahoz-Monfort, J., Pagel, J., Pescott, O. L., Schmucki, R., Simmonds, E. G. and O’Hara, R. B. (2020) Data Integration for Large-Scale Models of Species Distributions Trends in Ecology & Evolution, \emph{35}, 56--67

  Simpson, D. P., Rue, H., Riebler, A., Martins, T. G. and Sørbye, S. H. (2017) Penalising model component complexity: A principled, practical approach to constructing priors Statistical Science, \emph{32}, 1--28
}
