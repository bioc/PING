\name{RawReadsTrack}
\alias{RawReadsTrack}

\title{Raw reads density track}

\description{
This track displays the start position of the forward and reverse reads used as input for the segmentation step}

\usage{
RawReadsTrack(reads, chr, gen="gen")
}

\arguments{
\item{reads}{An object of class \code{\linkS4class{GRanges}}. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}

}

\value{
This function returns an object of class \code{DataTrack}.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{segmentPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{DataTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}