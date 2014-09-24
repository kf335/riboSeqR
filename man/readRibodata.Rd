\name{readRibodata}
\alias{readRibodata}
\title{
  Reads ribosomal and (optionally) rna data from alignment files.
}
\description{
  Reads flat text files (or compressed versions of these) containing
  strand, transcript name, start and sequence information for each alignment.
}
\usage{
readRibodata(riboFiles, rnaFiles, columns = c(strand = 1, seqname = 2,
start = 3, sequence = 4), zeroIndexed = TRUE, header = FALSE, replicates, seqnames)
}
\arguments{
  \item{riboFiles}{
    Filenames of ribosomal alignments.
}
  \item{rnaFiles}{
    Filenames of RNA alignments.
}
  \item{columns}{
    Columns of alignment files containing strand, transcript (seqname)
    name, start of alignment, and sequence.
  }
  \item{zeroIndexed}{
    Are the alignments zero-indexed (i.e., the first base in a sequence
    is 0). Defaults to TRUE, which will result in an adjustment to
    1-indexed data.}
  \item{header}{
    Does the alignment file have a header line? Defaults to FALSE.
}
  \item{replicates}{
    Replicate information for the files.
}
  \item{seqnames}{
    Transcript (seqname) names to be read into the object. 
}
}
\value{
\code{riboData} object.
}
\author{
Thomas J. Hardcastle
}

\examples{
datadir <- system.file("extdata", package = "riboSeqR")
ribofiles <- paste(datadir, 
                   "/chlamy236_plus_deNovo_plusOnly_Index", c(17,3,5,7), sep = "")
rnafiles <- paste(datadir, 
                  "/chlamy236_plus_deNovo_plusOnly_Index", c(10,12,14,16), sep = "")

riboDat <- readRibodata(ribofiles, rnafiles, replicates = c("WT", "WT", "M", "M"))
}
\keyword{IO}
\keyword{manip }
