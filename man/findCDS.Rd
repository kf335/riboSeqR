\name{findCDS}
\alias{findCDS}
\title{
Parses a transcriptome file looking for start/stop codons in frame.
}
\description{
Looks in the fastaFile for defined start and stop codons in frame with
one another. Reports the discovered coordinates as a GRanges object with
a `frame' value.
}
\usage{
findCDS(fastaFile, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))
}
\arguments{
  \item{fastaFile}{
    Fasta file of transcriptome sequences.
}
  \item{startCodon}{
    Vector of possible start codons. Defaults to "ATG". 
}
  \item{stopCodon}{
    Vector of possible stop codons. Defaults to c("TAG", "TAA", "TGA").
}
}

\value{
  A GRanges object.
}
\author{
  Thomas J. Hardcastle
}

\examples{
datadir <- system.file("extdata", package = "riboSeq")
chlamyFasta <- paste(datadir, "/rsem_chlamy236_deNovo.transcripts.fa", sep = "")
fastaCDS <- findCDS(fastaFile = chlamyFasta, 
                    startCodon = c("ATG"), 
                    stopCodon = c("TAG", "TAA", "TGA"))
}

\keyword{IO}
\keyword{manip}
