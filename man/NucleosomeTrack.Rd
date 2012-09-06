\name{NucleosomeTrack}
\alias{NucleosomeTrack}

\title{Track for the Nucleosome position}

\description{
This track displays the position of the predicted nucleosomes as well as the standard error associated with each prediction.
}

\usage{
NucleosomeTrack(PS, chr, gen="gen", filter=TRUE)
}

\arguments{
\item{PS}{A \code{data.frame} or \code{\linkS4class{pingList}}. The output of \code{PING} or \code{postPING}.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{filter}{A \code{logical}. Screen out abnormal predictions whose parameters are too far from others.}
}

\value{
This function returns an object of class \code{AnnotationTrack}.
}

\note{
\code{postPING} output is already filtered. Here, the filtering is done with default parameters. For a more customized filtering, run \code{postPING} with different values of sigmaB2.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{postPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{AnnotationTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}