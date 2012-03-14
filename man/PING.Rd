\name{ping}
\alias{ping}
\alias{PING}
\alias{ping-class}
\alias{code}
\alias{code,ping-method}
\alias{code,pingError-method}
\alias{code,pingList-method}
\alias{coerce}
\alias{coerce,pingList,RangedData-method}
\alias{coerce,pingList,data.frame-method}
\alias{coerce,RangedData,GenomeData-method}
\alias{coerce,data.frame,GenomeData-method}
\alias{coerce,data.frame,pingList-method}
\alias{coerce,data.frame,pingList-method}
\alias{density}
\alias{density,ping-method}
\alias{density,pingError-method}
\alias{score}
\alias{scoreForward}
\alias{scoreForward,ping-method}
\alias{scoreForward,pingList-method}
\alias{scoreForward,pingError-method}
\alias{scoreReverse}
\alias{scoreReverse,ping-method}
\alias{scoreReverse,pingList-method}
\alias{scoreReverse,pingError-method}
\alias{sigmaSqR}
\alias{sigmaSqR,ping-method}
\alias{sigmaSqR,pingList-method}
\alias{sigmaSqR,pingError-method}
\alias{sigmaSqF}
\alias{sigmaSqF,ping-method}
\alias{sigmaSqF,pingList-method}
\alias{sigmaSqF,pingError-method}
\alias{K}
\alias{K,ping-method}
\alias{K,pingList-method}
\alias{K,pingError-method}
\alias{plot}
\alias{plot,ping,segReads-method}
\alias{plot,pingError,segReads-method}
\alias{plot,pingList,pingList-method}
\alias{plot,data.frame,data.frame-method}
\alias{plot,pingList,segReadsList-method}
\alias{maxRange}
\alias{maxRange,ping-method}
\alias{maxRange,pingList-method}
\alias{maxRange,pingList-method}
\alias{maxRange,pingError-method}
\alias{minRange}
\alias{minRange,ping-method}
\alias{minRange,pingList-method}
\alias{minRange,pingList-method}
\alias{minRange,pingError-method}

% I think a few methods are missing, e.g. [[<- and [<-

\title{Estimation of binding site positions} 

\description{
This object contains Estimation of binding site positions and has the following  slots: 
segReadsList.
}

\usage{
PING(segReadsList)
}

\arguments{
\item{segReadsList}{This object contains segmentation of Genome}
}

\section{Methods}{
  \describe{
 \item{code}{\code{signature(x = "ping")}: return the error code for each list element (i.e. candidate region) of a PING object. If the string is empty, there were no errors.}
 \item{plot}{\code{signature(x = "ping")}: Plot all regions in the PING object. This might be long, and should only be used to plot a few regions, so subset the object before plotting.}
 \item{sigmaSqR}{\code{signature(x = "ping")}: return the variance parameter of the reverse (R) distribution for each binding event.}
 \item{sigmaSqF}{\code{signature(x = "ping")}: return the variance parameter of the forward (F) distribution for each binding event.}
 \item{score}{\code{signature(x = "ping")}: return the score for each binding event.} 
 \item{scoreF}{\code{signature(x = "ping")}: return the score  of the forward (F)  for each binding event.}
 \item{scoreR}{\code{signature(x = "ping")}: return the score  of the forward (R)  for each binding event.}
\item{maxRange}{\code{signature(x = "ping")}: return the range maximum.}
\item{minRange}{\code{signature(x = "ping")}: return the range minimal.}
 \item{K}{\code{signature(x = "ping")}: subset PING object.}
 \item{density}{\code{signature(x = "ping")}: return the density for each binding event.}
}
}

\author{
  Xuekui Zhang, \email{ubcxzhang@gmail.com} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Sangsoon Woo, \email{swoo@fhcrc.org}
}

\seealso{
\code{\link{ping}}
}
\keyword{models}