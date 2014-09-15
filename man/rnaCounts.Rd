\name{rnaCounts}
\alias{rnaCounts}
\title{
Extracts mRNA counts from a riboDat object for a set of coding sequence coordinates.
}
\description{
Takes mRNA count data from riboDat object, maps them to coding sequences
specified in GRanges object, and counts the total number of
hits. This is a crude approach intended to quickly produce comparable
data to ribosome footprint counts. More sophisticated alternatives,
addressing coverage variation, isoforms, multireads &c. have been widely
described in the literature on mRNA-seq analyses.
}
\usage{
rnaCounts(riboDat, CDS)
}

\arguments{
  \item{riboDat}{
    A \code{\link{riboData}} object containing the RNA-seq alignments.
}
  \item{CDS}{
    A \code{GRanges} object defining the coordinates of the coding
    sequences for which to acquire counts.
}
}
\details{
  The count data thus acquired can be compared to counts of ribosomal
  footprint data through a beta-binomial analysis (see vignette) to
  discover differential translation.
  }
\value{
  A matrix containing count data for the RNA-seq libraries.
}
\author{
Thomas J. Hardcastle
}

\seealso{
\code{\link{sliceCounts}}
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

# Extract counts of RNA hits from riboCount data.
rnaCounts <- rnaCounts(riboDat, ffCs@CDS)
}
\keyword{manip}
