
\name{summary.isdm}
\alias{summary.isdm}
\alias{print.summary.isdm}
\title{Summarizing integrated species distribution models (isdm)}
\usage{
  \method{summary}{isdm}(object, \dots)
  
  \method{print}{summary.isdm}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{object}{an object of class \code{"isdm"}, usually, a result of a
    call to \code{\link{isdm}}.}
  \item{x}{an object of class \code{"summary.isdm"}, usually, a result of a
    call to \code{summary.isdm}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{\dots}{not implemented.}
}
\description{
  These functions are all \code{\link{methods}} for class \code{isdm} or
  \code{summary.isdm} objects.
}
\details{
  The result summarises the posterior distribution of the model. For each component in turn, excluding those components not present in the current model. Some of the more usual descriptors are given, e.g. mean, standard deviation and some quantiles.
  
  Distribution descriptors are printed to within \code{"digits"} decimals places.
}
\value{
\code{summary.isdm} returns an object of class \code{"summary.isdm"}, a list with components
  
\item{DISTRIBUTION}{the posterior for the parameters in the distributional part of the \code{object}.}
\item{PO_BIAS}{the posterior for the bias terms (PO data) part of the \code{object}.}
\item{DC_ARTEFACT}{the posteriors for the artefacts relating to the DC component of \code{object}. Note that this may be NULL if no DC data were supplied. If not NULL, then this includes detectability parameters (those with \code{"alpha"} in the name).}
\item{AA_ARTEFACT}{the posteriors for the artefacts relating to the DC component of \code{object}. Note that this may be NULL if no AA data were supplied.}
\item{PA_ARTEFACT}{the posteriors for the artefacts relating to the PA component of \code{object}. Note that this may be NULL if no PA data were supplied.}
\item{SPATIAL}{the hyper-parameters for the spatial random effects included in \code{object}.}
\item{mark.lik}{the model evidence (marginal log-likelihood) of \code{object}.}
}
\seealso{
 \code{\link{isdm}}, \code{\link{predict.isdm}}.
}
