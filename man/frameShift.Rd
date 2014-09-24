\name{readingFrame}
\alias{readingFrame}
\alias{plotFS}
\title{
Analyses frame called ribosome footprint data within coding sequences
and identifies likely frame-shift of different length ribosome footprint reads.
}
\description{
Ribosome footprint data data can be used to identify the frame-shift
relative to start codon of the different n-mers. For each readlength
specified, the sum of alignments in the different frames is shown,
together with the maximum likelihood frame.
}
\usage{
readingFrame(coordinates, riboDat, rC, lengths = 26:30)
plotFS(fS, lengths, legend.text = c("Frame 0", "Frame 1", "Frame 2"), ...)
}
\arguments{
  \item{coordinates}{
    Coordinates of transcripts. Need not be given if a riboCoding object
    is specified in `rC'.
}
  \item{riboDat}{
    A \code{riboData} object containing the ribosome footprint
    data. Need not be given if a riboCoding object is specified in `rC'.
}
  \item{rC}{
    A riboCoding object.
}
  \item{lengths}{
    Lengths of reads to be analysed for frame-shift, or to be
    plotted. May be omitted in plotting, in which case all lengths will
    be plotted.
  }
  \item{fS}{
    The output of the readingFrame function, to be plotted.}
  \item{legend.text}{
    Text for legend.
  }
  \item{...}{
    Additional arguments to be passed to barplot function.
    }
}
\value{
  A matrix giving the number of aligned reads in each frame for each
  length, and the maximum likelihood frame.
}

\author{
Thomas J. Hardcastle
}

\seealso{
\code{\link{frameCounting}}
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
rCs <- frameCounting(riboDat, fastaCDS)

# analysis of frame shift for 27 and 28-mers.
fS <- readingFrame(rC = rCs, lengths = 27:28)
plotFS(fS)

}
\keyword{manip}
