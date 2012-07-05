\name{setParaPriorPING}
\alias{setParaPriorPING}

\title{
	Function that returns a list of parameters that can be used as an argument of PING or postPING.
}

\description{
	This function takes from 0 to 6 parameters as argument, check if they are valid and returns a list to be used in a call to PING or postPING. The default parameters are set for MNase data.
}

\usage{
	setParaPriorPING(xi=150,rho=0.8,alpha=20,beta=20000,lambda=-0.000064,dMu=200)
}

\arguments{
	\item{xi}{ An \code{integer}. The average DNA fragment size.}
	\item{rho}{ An \code{integer}. A variance parameter for the average DNA fragment size distribution.}
	\item{alpha}{ An \code{integer}. First hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model}
	\item{beta}{ An \code{integer}. Second hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model}
	\item{lambda}{ An \code{integer}. The precision of the prior for mu used for histone data.}
	\item{dMu}{ An \code{integer}. Our best guess for the distance between two neighboring nucleosomes.}
}

\value{
	Returns a list of 6 parameters to be used in PING.
}

\author{
Renan Sauteraud
}

\seealso{
	\code{PING}
	\code{postPING}
}

\examples{
	# set prior for MNase data
	paraPrior<-setParaPriorPING(xi=150,rho=0.8,alpha=20,beta=20000,lambda=-0.000064,dMu=200) 
	# set prior for sonication data
	paraPrior<-setParaPriorPING(xi=150,rho=1.2,alpha=10,beta=20000,lambda=-0.000064,dMu=200) 
}

\keyword{functions}
