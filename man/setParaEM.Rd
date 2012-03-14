\name{setParaEM} 
\alias{setParaEM} 
\alias{paraEMH} 
\alias{paraEMTF} 
\title{Set convergence parameters of the EM algorithm} 

\description{ 
This function can be used to change the internal PING parameters for the EM algorithm. This function should only be called if you really now what you are doing!.}

\usage{
setParaEM(minK=0,maxK=0,tol=1e-4,B=100,mSelect="AIC3",mergePeaks=TRUE,mapCorrect=TRUE)
}

\arguments{
  \item{minK}{An integer value. The minimum number of binding events per region. If value is 0, the minimum number is automatically calculated.}
  \item{maxK}{An integer value. The maximum number of binding events per region. If value is 0, the maximum number is automatically calculated.}
  \item{tol}{The tolerance for the EM algorithm}
  \item{B}{An integer value. The maximum number of iterations to be used.}
  \item{mSelect}{A character string specifying the information criteria to be used when selecting the number of binding events.}
  \item{mergePeaks}{A logical value stating whether overlapping binding events should be picked.}
  \item{mapCorrect}{Should mappability profiles be incorporated in the estimation, that is missing reads estimated.}
}


\value{
No value returned. The function simply modifies the internal variables 'paraEMH'.
}
\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" GenomeBiology, under review.
}
\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\seealso{setParaPrior}
\examples{
# Using mSelect="BIC"
setParaEM(minK=1,maxK=8,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE)
# Using mSelect="AIC"
setParaEM(minK=1,maxK=8,tol=1e-4,B=100,mSelect="AIC",mergePeaks=TRUE,mapCorrect=TRUE)
}
\keyword{data}