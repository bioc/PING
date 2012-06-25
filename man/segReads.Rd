\name{segReads}
\alias{segReads}
\alias{segReads-class}
\alias{map}
\alias{map,segReads-method}
\alias{map,segReadsList-method}
\alias{\%in\%,ANY,segReads-method}


\title{Segment the genome into candidate regions} 

\description{
Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PING.
}

\usage{
segReads(yF, yR, cF, cR, map, chr)
}

\arguments{
\item{yF}{This object contains an ExpressionSet}
\item{yR}{ String containing the genome name used (vector).}
\item{cF}{ String containing the name of chromosome used (vector).}
\item{cR}{ String containing the Position of the sequences (vector).}
\item{map}{ String containing the copy number of sequence (vector).}
\item{chr}{ String containing the expresion data of enriched region (matrix with n column.}

}

\section{Methods}{
  \describe{
 \item{map}{\code{signature(x = "ping")}: subset PING object.}
}
}


\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" PlosONE, under review.
}
\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\seealso{
\code{\link{ping}}
}

\keyword{models}