\name{filterHits}
\alias{filterHits}
\title{
  Filters framecalled data based on the mean number of hits observed
  across a replicate group.
}
\description{
If no ribosomal footprints of the correct lengths (and frames) are seen
at a coding sequence in any replicate group, this sequence is unlikely
to be translated (and is therefore likely to be uninteresting). This
function filters out these coding sequences, leaving only those with a
minimum number of hits in at least one replicate group, and a minimum
number of unique sequences aligning in at least one replicate group (to
exclude single stacks of sequenced reads passing the filtering).
}
\usage{
filterHits(fCs, lengths = 27, frames, hitMean = 10, unqhitMean = 1)
}
\arguments{
  \item{fCs}{
    \code{riboCoding} object defining the number of ribosome footprint
    reads observed over a set of coordinates, their lengths, and their
    frame relative to coding start sites.
}
  \item{lengths}{
    The lengths of ribosome footprint reads to be used in filtering.
}
  \item{frames}{
    If given, the frames of the ribosome footprint reads to be used in
    filtering. Should be of equal length to the `lengths' parameter -
    see examples.
}
  \item{hitMean}{
    The mean number of hits within a replicate group that should be
    observed to pass filtering.
}
  \item{unqhitMean}{
    The mean number of unique sequences within a replicate group that
    should be observed to pass filtering. This parameter is intended to
    avoid cases where a coding sequence is deemed to be expressed based
    on a few highly expressed sequences.
}
}
\details{
  Frames can be given as a single vector (which specifies the frames
  used for all lengths of footprints, or as a list, specifying the frame
  for each length given in `lengths'.  
}
\value{
  \code{riboCoding} object containing the filtered coding sequences and
  the associated numbers of ribosome footprint reads.
}
\author{
  Thomas J. Hardcastle
}

\seealso{
  \code{frameCounting}
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

# filter coding sequences. 27-mers are principally in the 0-frame,
# 28-mers are principally in the 2-frame relative to coding start (see
# readingFrame function).

ffCs <- filterHits(fCs, lengths = c(27, 28), frames = list(0, 2), 
                   hitMean = 50, unqhitMean = 10)


}
\keyword{manip}
