.prepareSlice <-
  function(slices, array) {
    slice <- paste(
               paste(
                 sapply(slices, function(x)
                        if(!is.null(x)) {
                          wdif <- which(diff(x) != 1)
                          paste("c(", paste(sapply(lapply(split(x, rep(1:(length(wdif) + 1), c(wdif, length(x)) - c(1, wdif + 1) + 1)), range), paste, collapse = ":"), collapse = ","), ")", sep = "")
                        } else return("")),                      
                 collapse = ","),
               paste(rep(",", length(dim(array)) - length(slices)), collapse = ""), sep = "")
    slice
  }

.sliceArray <-
  function(slices, array, drop = FALSE) {
    if((is.vector(array) & sum(!sapply(slices, is.null)) > 1) || (is.array(array) & length(slices) > length(dim(array)))) warning("dimensions of slice exceed dimensions of array")
    if(is.vector(array)) return(array[slices[[min(which(!sapply(slices, is.null)))]]])
    dropText <- c(", drop = FALSE", ", drop = TRUE")[as.numeric(drop) + 1]
    slice <- .prepareSlice(slices, array)
    sarray <- eval(parse(text = paste("array[", slice, dropText, "]", sep = "")))
    sarray
  }
