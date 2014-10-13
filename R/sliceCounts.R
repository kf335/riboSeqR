sliceCounts <-
function(rC, lengths = 27, frames) {  
  if(missing(frames)) frames <- list(0:2)
  if(!is.list(frames)) frames <- list(frames)
  if(length(frames) == 1) frames <- rep(frames, length(lengths))
  if(length(lengths) == 1) lengths <- rep(length, length(frames))
  if(length(frames) != length(lengths)) stop("Length of 'lengths' parameter and 'frames' parameter must be equal (or one) if given")

  fH <- do.call("cbind", lapply(1:ncol(rC@hits), function(jj) {
    (do.call("+", lapply(1:length(lengths), function(ll) {
      rC@hits[,jj,as.character(frames[[ll]]), as.character(lengths[ll])]
    })))
  }))

  
  fH

}

rnaCounts <- function(riboDat, CDS) {
  z <- do.call("cbind", lapply(riboDat@rnaGR, function(x) table(subjectHits(findOverlaps(x, CDS)))))
}
