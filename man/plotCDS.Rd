\name{plotCDS}
\alias{plotCDS}

\title{
Plots average ribosome footprint alignment to coding sequences at 5' and
3' ends.
}
\description{
For each sample, the average (normalised by translation abundance over
transcript) of the ribosome footprints of a given length alignments at
the 5' and 3' ends of all specified transcripts beginning at each base
relative to coding start/end are plotted. The bases are colour coded
relative to start codon.
}
\usage{
plotCDS(coordinates, riboDat, lengths = 27, min5p = -20, max5p = 200,
min3p = -200, max3p = 20, cap, main = "", plot = TRUE, ...)
}
\arguments{
  \item{coordinates}{
    Coordinates (as a \code{GRanges} object) of the coding sequences.
}
  \item{riboDat}{
    \code{\link{riboData}} object containing ribosome footprint data.
}
  \item{lengths}{
    Lengths of footprints to be plotted. May be given as a vector, in
    which case multiple plots will be produced.
}
  \item{min5p}{
    The distance upstream of the translation start to be plotted.
}
  \item{max5p}{
    The distance downstream of the translation start to be plotted.
}
  \item{min3p}{
    The distance upstream of the translation end to be plotted.
}
  \item{max3p}{
    The distance downstream of the translation end to be plotted.
  }
  \item{cap}{
    If given, caps the height of plotted values.
    }
  \item{main}{
    Title of the plot.
}
  \item{plot}{
    Should the acquired matrix of mean expression be plotted? Defaults
    to TRUE.
  }
  \item{...}{
    Additional arguments to be passed to 'plot' and 'axes'.
  }
}

\value{
Invisibly returns lists of lists of matrices containing weighted
averages plotted for each sample/length combination.
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

# CDS coordinates
chlamyFasta <- paste(datadir, "/rsem_chlamy236_deNovo.transcripts.fa", sep = "")
fastaCDS <- findCDS(fastaFile = chlamyFasta, 
                    startCodon = c("ATG"), 
                    stopCodon = c("TAG", "TAA", "TGA"))

# frame calling
fCs <- frameCounting(riboDat, fastaCDS)

# analysis of frame shift for 27 and 28-mers.
fS <- readingFrame(rC = fCs, lengths = 27:28)

# filter coding sequences. 27-mers are principally in the 0-frame,
# 28-mers are principally in the 2-frame relative to coding start (see
# readingFrame function).

ffCs <- filterHits(fCs, lengths = c(27, 28), frames = list(0, 2), 
                   hitMean = 50, unqhitMean = 10, fS = fS)

plotCDS(coordinates = ffCs@CDS, riboDat = riboDat, lengths = 27)

}
\keyword{hplot}
