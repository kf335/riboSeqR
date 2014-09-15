\name{sliceCounts}
\alias{sliceCounts}
\title{
Slices out count data from riboCoding object for use in differential
translation analyses.
}
\description{
For any given coding sequence, multiple lengths of reads in various
frames (relative to coding start) may align. This function extracts
specific size-classes and frames of ribosome footprint reads and sums
them to give a single value for each coding sequence and each sequencing
library, for use in downstream analysis.
}
\usage{
sliceCounts(rC, lengths = 27, frames)
}
\arguments{
  \item{rC}{
    A \code{\link{riboCoding}} object containing the coordinates of the
    coding sequences and the number of ribosomal footprint reads of
    various classes to be found in each.
}
  \item{lengths}{
    Lengths of ribosome footprints to inform count data.
}
  \item{frames}{
    Frames of ribosome footprints (relative to coding start site). If
    omitted, all frames are used. 
}
}
\details{
  Frames can be given as a single vector (which specifies the frames
  used for all lengths of footprints, or as a list, specifying the frame
  for each length given in `lengths'.

  The count data thus acquired can be compared to counts of RNA-seq data
  through a beta-binomial analysis (see vignette) to discover
  differential translation.
}
\value{
A matrix containing counts of ribosomal footprint matches to coding
sequences specified in riboCoding object `rC'.
}
\author{
Thomas J. Hardcastle
}

\seealso{
\code{\link{rnaCounts}}
}
\examples{
#ribosomal footprint data
datadir <- system.file("extdata", package = "riboSeq")
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

# Extract counts of ribosomal footprints from riboCount data.

riboCounts <- sliceCounts(ffCs, lengths = c(27, 28), frames = list(0,
2))

}
\keyword{manip}
