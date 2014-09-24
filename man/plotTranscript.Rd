\name{plotTranscript}
\alias{plotTranscript}
\title{
Plots ribosome footprint abundance and mRNA coverage (if available) for
a specific transcript.
}
\description{
Abundances of ribosomal footprints of a given size class are plotted on
a transcript. The footprints are colour coded according to the first
base of the transcript, and not any coding start site, to allow for
multiple coding start sites on a given transcript. Coding regions may
simultaneously be plotted and colour coded under the same scheme.
}
\usage{
plotTranscript(transcript, coordinates, annotation, riboData, length =
27, frameShift = 0, cap, riboScale, rnaScale, baseLim, main, note = "", ...)
}
\arguments{
  \item{transcript}{
    The name of the transcript to be plotted.
  }
  \item{coordinates}{
    A GRanges object containing any coding regions on the transcript.
  }
  \item{annotation}{
    A GRanges object containing annotated coding coordinates to be
    plotted as bars above the figure.
  }
  \item{riboData}{
    A \code{\link{riboData}} object containing the ribosome footprint
    (and optionally, RNA-seq) data.
  }
  \item{length}{
    Size class of ribosome footprint data to be plotted.
  }
  \item{frameShift}{
    Frameshift for the ribosome footprint data. See Details.
  }
  \item{cap}{
   Cap on the largest value that will be plotted as an abundance of the
   ribosome footprint data.
  }
  \item{riboScale}{
    Scale to be used on the ribosome footprint axis.
  }
  \item{rnaScale}{
    Scale to be used on the RNA-seq coverage axis.
  }
  \item{baseLim}{
    Limits of the bases of the transcript to be plotted (i.e., the
    x-axis). If missing, the full transcript will be plotted.
  }
  \item{main}{
    Optional title for the plot.
    }
  \item{note}{
    Additional note to be added to plot titles (in addition to transcript
    and sample names).
  }
  \item{...}{
    Additional arguments to be passed to plotting function.
    }
}

\details{
  The readingFrame value allows the colour-coding of the ribosome
  footprints to be shifted so that the colours of the coding sequences match
  the colours of the ribosome footprint data. E.g., if 28-mers are
  predominantly in frame 2 relative to coding start, a value of
  `readingFrame=2' will ensure that 28-mers in a coding region will take
  the same colour as that coding region if they are in the correct
  relative frame.
}
\value{
NULL; plotting function.
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

plotTranscript("CUFF.37930.1", coordinates = ffCs@CDS, 
               riboData = riboDat, length = 27, cap = 200)               

}
\keyword{hplot}
