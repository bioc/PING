\name{NucScoreTrack}
\alias{NucScoreTrack}

\title{Track for the Nucleosome score}

\description{
This track displays an histogram of the score associated with the predicted nucleosomes.
}

\usage{
NucScoreTrack(PS, chr, gen="gen", ...)
}

\arguments{
\item{PS}{A \code{data.frame} or \code{\linkS4class{pingList}}. The output of \code{ipostPING} or \code{PING}.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{...}{A list of supplemental argument that can be provided to the \code{\link[Gviz]{DataTrack}} track constructor.}
}

\value{
This function returns an object of class \code{DataTrack}.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{postPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{DataTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}
