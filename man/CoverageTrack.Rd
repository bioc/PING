\name{CoverageTrack}
\alias{CoverageTrack}

\title{Reads coverage track}

\description{
This track displays the coverage of the genomic region by the reads used as input for the segmentation step.}

\usage{
CoverageTrack(reads, chr, gen="gen", FragmentLength=200)
}

\arguments{
\item{reads}{An object of class \code{\linkS4class{GRanges}}. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{FragmentLength}{An \code{integer}. The fixed length of the reads.}
}

\value{
This function returns an object of class \code{AnnotationTrack}.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{postPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{AnnotationTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}