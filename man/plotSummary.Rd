\name{plotSummary}
\alias{plotSummary}
\alias{CoverageTrack}
\alias{RawReadsTrack}
\alias{NucleosomeTrack}
%Ultimately, I'll make the doc for each track

\title{
Plot a summary of the prediction for given ranges.
}

\description{
This function use Gviz package to summarize the nucleosome position prediction from postPING.
It displays 4 tracks:
A \code{CoverageTrack} which displays the coverage of the given region by the reads used as input of the segmentation step.
A \code{RawReadsTrack} which displays the start position of the forward and reverse reads.
A \code{NucleosomeTrack} which displays the position of the predicted nucleosomes as well as the standard error associated with the prediction.
}

\usage{
plotSummary(PS, reads, chr, gen="gen", from=NULL, to=NULL, FragmentLength=200, title="")
}

\arguments{
\item{PS}{A \code{data.frame}. The output of the \code{postPING}.}
\item{reads}{A \code{GRanges} object. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{A character string for the genome name.}
\item{from}{An \code{integer}, the first base of the plot}
\item{to}{An \code{integer}, the last base of the plot}
\item{FragmentLength}{An \code{integer}, the length of XSET profile extension.}
\item{title}{An optional character string that can be prepend to the automatically generated title.}
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

\keyword{functions}