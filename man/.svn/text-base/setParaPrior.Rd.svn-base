\name{setParaPrior} 
\alias{setParaPrior} 
\alias{paraPriorTF}
\alias{paraPriorH}

\title{Set convergence parameters of the EM algorithm} 

\description{ 
This function can be used to change the internal PING parameters for the prior distribution. This function should only be called if you really now what you are doing! In particular, you may want to specify the average DNA fragment size for your sample by changing the 'xi' parameter.
}
\usage{
setParaPrior(xi=150,rho=1.2,alpha=10,beta=20000,lambda=-0.000064,dMu=200) 

}
\arguments{
  \item{xi}{Our best guest for the average DNA fragment size.}
  \item{rho}{A variance parameter for the average DNA fragment size distribution.}
  \item{alpha}{First hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model.}
  \item{beta}{First hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model.}
  \item{lambda}{The precision of the prior for mu used for histone data.}
  \item{dMu}{Our best guess for the distance between two neighboring nucleosomes.}
}
%\details{
%}
\value{
No value returned. The function simply modifies the internal variables 'paraPriorH'.
}
\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" GenomeBiology, under review.
}
\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\seealso{setParaEM}
\examples{
# set prior for MNase data
setParaPrior(xi=150,rho=0.8,alpha=20,beta=20000,lambda=-0.000064,dMu=200) 
# set prior for sonication data
setParaPrior(xi=150,rho=1.2,alpha=10,beta=20000,lambda=-0.000064,dMu=200) 
}

\keyword{data}