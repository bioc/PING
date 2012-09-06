\name{plotSummary}
\alias{plotSummary}

\title{
Plot a summary of the prediction for given ranges.
}

\description{
This function use Gviz package to summarize the nucleosome position prediction from postPING.
}

\usage{
plotSummary(PS, reads, chr, gen="gen", from=NULL, to=NULL, FragmentLength=200, title="", filter=TRUE)
}

\arguments{
\item{PS}{An object of class \code{\linkS4class{pingList}} or \code{data.frame}. The output of \code{PING} or \code{postPING}. If a list of \code{data.frame} or \code{\linkS4class{pingList}} is passed, multiple \code{NucleosomeTrack} will be plotted.}
\item{reads}{A \code{GRanges} object. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{from}{An \code{integer}, the first base of the plot}
\item{to}{An \code{integer}, the last base of the plot}
\item{FragmentLength}{An \code{integer}, the length of XSET profile extension in the \code{CoverageTrack}.}
\item{title}{An optional character string that can be prepend to the automatically generated title.}
\item{filter}{A \code{logical}. Screen out abnormal predictions whose parameters are too far from others in the \code{NucleosomeTrack}.}
}

\note{
This function cannot be used for paired-end sequencing data as the reads have to be in GRanges format.
}

\value{
This function returns an invisible \code{list} of tracks. The tracks can be used in other Gviz plots.
}

\author{
Renan Sauteraud
}

\seealso{
\code{\link{NucleosomeTrack}}
\code{\link{RawReadsTrack}}
\code{\link{CoverageTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}