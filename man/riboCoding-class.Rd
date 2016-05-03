\name{riboCoding-class}
\docType{class}
\alias{riboCoding}
\alias{riboCoding-class}
\alias{[,riboCoding-method}
\alias{[,riboCoding,ANY-method}
\alias{show,riboCoding-method}

\title{Class "riboCoding"}

\description{
  The \code{riboCoding} class contains a set of coordinates defining
  coding sequences, a set of replicate data for the experimental
  samples, an array of ribosome footprint abundances for each
  coding sequence split by size class and frame, and a similar array
  describing the abundance of unique sequences aligning within each
  coding sequence. 
}

\section{Slots}{
  \describe{
    \item{\code{CDS}:}{Object of class \code{"GRanges"} defining the
      coordinates of coding sequences.}
    \item{\code{hits}:}{An array describing the abundances of ribosome
      footprints, split by size class and frame (relative to coding start)
      for each of the coding sequences defined in the `CDS' slot.}
    \item{\code{unqHits}:}{An array describing the abundances of unique
      sequences of ribosome footprints, split by size class and frame
      (relative to coding start) for each of the coding sequences defined in
      the `CDS' slot.}
    \item{\code{replicates}:}{A factor defining the replicate structure
      of the samples for which data are held. Replicate samples will have
      the same level in this factor.}
    }
}

\section{Methods}{
  Methods `[' and `show' are defined for this class.
}

\author{Thomas J. Hardcastle}

\keyword{classes}
