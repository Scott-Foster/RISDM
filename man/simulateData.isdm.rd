\name{simulateData.isdm}
\alias{simulateData.isdm}
\title{Simulates data of different types from an underlying point process.}
\description{This function simulates some data based on a number of fairly strict assumptions. }
\usage{
 simulateData.isdm( expected.pop.size=400000, expected.n.PO=300, n.PA=150, n.AA=50, n.DC=50,
                           coefs=list(dist=c(-1,-0.25,0.75), bias=c(-2,-1.5)), 
                           DC.pis=matrix( c(0.8,0.76, 0.7,0.73, 0.82,0.67), nrow=3, ncol=2, byrow=TRUE),
                           transect.size = 0.125, #a proportion of cell size.
                           rasterBoundary=NULL,
                           control=list())
}
\arguments{
\item{expected.pop.size}{coming...}
\item{expected.n.PO}{coming...}
\item{n.PA}{coming...}
\item{n.AA}{coming...}
\item{n.DC}{coming...}
\item{coefs}{coming...}
\item{bias}{coming...}
\item{DC.pis}{coming...}
\item{transect.size}{coming...}
\item{rasterBoundary}{coming...}
\item{control}{coming...}
}
\details{
 This function generates some fake data. It was written largely for internal testing purposes, but is made available just in case others find it useful (or useful to hack).
}
\value{
 NULL
}

\seealso{\code{\link{isdm}}, \code{\link{predict.isdm}}}

\author{Scott D. Foster}

