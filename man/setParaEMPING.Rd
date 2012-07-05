\name{setParaEMPING}
\alias{setParaEMPING}

\title{
	Function that returns a list of parameters for the EM algorithm that can be used as an argument of PING.
}

\description{
	This function takes from 0 to 7 EM algorithm parameters as argument, check if they are valid and returns a list to be used in a call to PING.
}

\usage{
	 setParaEM(minK=0,maxK=0,tol=1e-4,B=100,mSelect="AIC3",mergePeaks=TRUE,mapCorrect=TRUE)
}

\arguments{
	\item{minK}{ An \code{integer}. The minimum number of binding events per region. If the value is 0, the minimum number is automatically calculated.}
	\item{maxK}{ An \code{integer}. The maximum number of binding events per region. If the value is 0, the maximum number is automatically calculated.}
	\item{tol}{ A \code{numeric}. The tolerance for the EM algorithm.}
	\item{B}{ An \code{integer}. The maximum number of iterations to be used.}
	\item{mSelect}{ A character string specifying the information criteria to be used when selecting the number of binding events.}
	\item{mergePeaks}{ A \code{logical} stating whether overlapping binding events should be picked.}
	\item{mapCorrect}{ A \code{logical} stating whether mappability profiles should be incorporated in the estimation, i.e: missing reads estimated.}
}

\value{
	Returns a list of 7 parameters to be used in PING.
}

\author{
Renan Sauteraud
}

\seealso{
	\code{PING}
	\code{postPING}
}

\keyword{functions}
