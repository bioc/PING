\name{bam2gr}
\alias{bam2gr}

\title{pre-process bam files}

\description{
Reads a bam file using \code{Rsamtools} and extract the reads for each chromosome.
}

\usage{
bam2gr(bamFile, chr=NULL)
}

\arguments{
\item{bamFile}{A character string, the name of the .bam file to read.}
\item{chr}{An optional character string. If specified, only the selected chromosome will be returned. Speed up the computation.}
%\item{outpath}{A character string, the path of the folder where the files should be created if save=TRUE}
%\item{save}{A \code{logical}. If set to TRUE, the function creates a file for each chromosome found in \code{bamFile}. If set to FALSE, returns a list of lists that can be used as input for the \code{segmentPING} function.}
}

\value{
Returns a \code{GRanges} of all the reads for each chromosome. 
}

\author{Renan Sauteraud}

\seealso{
\code{\link{segmentPING}}
}
\keyword{functions}
