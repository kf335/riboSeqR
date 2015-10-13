\name{lengthDist}
\alias{lengthDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function plotting the density of lengths in a riboData object.
}
\description{
riboData object contain data on multiple lengths of sequenced
reads. This function plots the density of those lengths.
}
\usage{
lengthDist(riboDat, add = FALSE, legend = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{riboDat}{
    riboData object to be plotted
}
  \item{add}{
    Should the curve be added to an existing plot, or a new plot drawn?
}
  \item{legend}{
    Legend text, if given.
}
  \item{\dots}{
    Additional arguments to be passed to `plot' or `legend'.
}
}
\author{
Thomas J. Hardcastle
}

\examples{
#ribosomal footprint data
datadir <- system.file("extdata", package = "riboSeqR")
ribofiles <- paste(datadir, 
                   "/chlamy236_plus_deNovo_plusOnly_Index", c(17,3,5,7), sep = "")
rnafiles <- paste(datadir, 
                  "/chlamy236_plus_deNovo_plusOnly_Index", c(10,12,14,16), sep = "")

riboDat <- readRibodata(ribofiles, rnafiles, replicates = c("WT", "WT",
"M", "M"))

lengthDist(riboDat)
}

\keyword{plot}
