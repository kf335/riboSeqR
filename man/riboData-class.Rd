\name{riboData-class}
\docType{class}
\alias{riboData}
\alias{riboData-class}

\title{Class "riboData"}

\description{
  The \code{riboData} class contains a list of GRanges objects
  containing ribosome footprint alignment data, a factor defining the
  replicate structure of the samples involved, and (optionally) a list
  of GRanges objects containing RNA-seq alignment data (paired with the
  ribosome footprint data). It will generally be created by the
  `readRibodata' function and not directly by the user.
}

\section{Slots}{
  \describe{
    \item{\code{riboGR}:}{List of \code{"GRanges"} objects (one for each
      sequenced sample) describing the alignments of ribosomal footprint
      data to a transcriptome.}
    \item{\code{rnaGR}:}{List of \code{"GRanges"} objects (one for each
      sequenced sample, paired with `riboGR' slot) describing the
      alignments of RNA-seq data to a transcriptome.}
    \item{\code{replicates}:}{A factor defining the replicate structure
      of the samples for which data are held. Replicate samples will have
      the same level in this factor.}
    }
}

\section{Methods}{
  No methods are currently defined for this class.
}

\author{Thomas J. Hardcastle}

\keyword{classes}
