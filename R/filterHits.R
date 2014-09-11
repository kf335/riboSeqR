filterHits <-
function(fCs, lengths =  27, frames, hitMean = 10, unqhitMean = 1) {
  if(missing(frames)) frames <- list(0:2)
  if(!is.list(frames)) frames <- list(frames)
  if(length(frames) == 1) frames <- rep(frames, length(lengths))
  if(length(lengths) == 1) lengths <- rep(lengths, length(frames))
  if(length(frames) != length(lengths)) stop("Length of 'lengths' parameter and 'frames' parameter must be equal (or one) if given")

  fH <- do.call("cbind", lapply(levels(fCs@replicates), function(rep) {
    rowMeans(Reduce("+", lapply(1:length(lengths), function(ll) {
        fCs@hits[,which(fCs@replicates == rep),as.character(frames[[ll]]), as.character(lengths[ll]), drop = FALSE]
      }))) > hitMean
  }))

  ufH <- do.call("cbind", lapply(levels(fCs@replicates), function(rep) {
    rowMeans(Reduce("+", lapply(1:length(lengths), function(ll) {
      fCs@unqHits[,which(fCs@replicates == rep),as.character(frames[[ll]]), as.character(lengths[ll]), drop = FALSE]
    }))) > unqhitMean
  }))

  filtration <- which(rowSums(fH, na.rm = TRUE) > 0 & rowSums(ufH, na.rm = TRUE) > 0)
  fCs[filtration,]
}
