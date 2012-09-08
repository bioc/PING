\name{prePING}
\alias{prePING}

\title{pre-process bam files}

\description{
Reads a bam file using \code{Rsamtools} and extract the reads for each chromosome.
}

\usage{
prePING(bamFile, outpath="./", save=TRUE)
}

\arguments{
\item{bamFile}{A character string, the name of the .bam file to read.}
\item{outpath}{A character string, the path of the folder where the files should be created if save=TRUE}
\item{save}{A \code{logical}. If set to TRUE, the function creates a file for each chromosome found in \code{bamFile}. If set to FALSE, returns a list of lists that can be used as input for the \code{segmentPING} function.}
}

\value{
Returns a \code{list} with one attribute for each chromosome. 
}

\author{Renan Sauteraud}

\seealso{
\code{\link{segmentPING}}
}
\keyword{functions}