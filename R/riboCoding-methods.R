setMethod("show", "riboCoding", function(object) {
  show(object@CDS)

  cat('\nSlot "replicates"\n')
  print(object@replicates)
  
  cat('\nHits array; dimension ', dim(object@hits), '\n')
  cat('\nUnique hits array; dimension ', dim(object@unqHits), '\n')  
})

setMethod("[", signature(x = "riboCoding", i = "ANY", j = "ANY"), function(x, i, j, ..., drop = FALSE) {
#  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(!missing(i))
    {
      i <- as.vector(i)
      if(is.logical(i)) i <- which(i)
      x@CDS <- x@CDS[i]      
    }  
  if(!missing(j))
    {
      j <- as.vector(j)
      if(is.logical(j)) j <- which(j)
      x@replicates <- x@replicates[j]
    }
  if(missing(i)) i <- NULL
  if(missing(j)) j <- NULL
  
  x@hits <- .sliceArray(list(i, j), x@hits)
  x@unqHits <- .sliceArray(list(i, j), x@unqHits)
  x
})
