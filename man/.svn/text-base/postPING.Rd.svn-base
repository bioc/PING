\name{postPING} 
\alias{postPING} 


\title{Post process Estimation of binding site positions obtained from PING} 
\description{ 
   Post process Estimation of binding site positions obtained from PING. Refit mixture models with stronger prior in candidate regions contain potential problems, and then convert final result into dataframe.}
\usage{
postPING(ping, seg, rho=8, sigmaB2=2500,rho1=0.8,alpha1=20,alpha2=100,beta2=100000,xi=150, min.dist=100, lambda=-0.000064)
}


\arguments{
  \item{ping}{A 'pingList' object containing estimation of nuclesome positions, result of 'PING' function.}
  \item{seg}{An object of class 'segmentReadsList' containing the results for all regions pre-processed, 'segmentReads' function.}
  \item{rho1, sigmaB2, rho, alpha1, alpha2, beta2, xi}{Integer values, the parameters in the prior of mixture models to be re-fitted.}
  \item{min.dist}{The minimum distance of two adjacent nucleosomes predicted from different candidate regions, smaller than that will be treated as duplicated predictions for the same nucleosomes.}
  \item{lambda}{The lambda control Gaussian Markov Random Field prior on the distance of adjacent nucleosomes, we do not recommend user change the default value.}
}

\value{
An dataframe containing the estimation of binding site positions. 
}
\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" PlosONE, under review.
}
\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\note{ Based on our experiemt on a few real data sets, we suggestion to use following values of parameters. For sonication data we use {rho1=1.2; sigmaB2=6400;rho=15;alpha1=10; alpha2=98; beta2=200000}.
For MNase data we use {rho1=3; sigmaB2=4900; rho=8; alpha1=20; alpha2=100; beta2=100000}. The value of xi depends on specy of sample, since that affect the length of linker-DNA. For example, we use xi=160 for yeast and xi=200 for mouse.
}
\keyword{data}
\keyword{models}