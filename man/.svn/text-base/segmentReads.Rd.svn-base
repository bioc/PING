\name{segmentReads} 
\alias{segmentReads} 


\title{Segment the genome into candidate regions} 
\description{ 
    Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PING.}
\usage{
segmentReads(data, dataC=NULL, map=NULL,minReads=2, minReadsInRegion=3,jitter=FALSE, maxLregion=1200,minLregion=80)
}


\arguments{
  \item{data}{A 'GenomeData' object containing the IP reads. See details for more information on how to set up the data.}
  \item{dataC}{A 'GenomeData' object containing the control reads. Set to NULL by default, i.e. no control.}
  \item{map}{A 'RangedData' object containing the mappability profiles. Set to NULL by default, i.e. no profiles.}
  \item{minReads}{The minimum number of F/R reads to be present in the sliding window. Set to NULL by default, and this number is automatically calculated.}
  \item{minReadsInRegion}{The minimum number of F/R reads to be present in the region.}
  \item{jitter}{A logical value stating whether some noise should be added to the read locations. This is recommended if the read positions have lots of duplicates.}
\item{maxLregion}{The maximum length. We do not suggest to use too long region, since that need to fit mixture models with too many components.}
\item{minLregion}{The minimum length.}
}

\value{
An object of class 'segmentReadsList' containing the results for all regions pre-processed. 
}
\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" PlosONE, under review.
}
\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\seealso{segmentReads}
\examples{
# Read data
path<- system.file("extdata",package="PING")
## Note that the col name for the chromosome needs to be space and not chr
dataIP<-read.table(file.path(path, "GSM351492_R4_chr1.bed"),header=TRUE,colClasses=c("factor","integer","integer","factor"))
dataIP<-as(dataIP,"RangedData")
dataIP<-as(dataIP,"GenomeData")

seg<-segmentReads(dataIP, minReads=2)
}

\keyword{data}
\keyword{models}
