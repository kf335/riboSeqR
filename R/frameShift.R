frameShift <- function(...) {
  stop("This function is obsolete. Please use readingFrame instead.")
}

readingFrame <- function(coordinates, riboDat, rC, lengths = 26:30) {    
  if(missing(coordinates)) coordinates <- NULL
  if(missing(riboDat)) riboDat <- NULL
  if(missing(rC)) rC <- NULL
  
  frameCounts <- sapply(lengths, function(length) {
    if((is.null(coordinates) | is.null(riboDat)) & !is.null(rC)) {
      frameCounts <- rowSums(
                       sapply(0:2, function(frame) {
      frCDS <- which(rC@CDS$frame == frame)
      if(length(frCDS) == 0) return(c(0, 0, 0))
      z <- .sliceArray(list(frCDS), rC@hits, drop = FALSE)
      zz <- z[,,,as.character(length),drop= FALSE]
      cs <- colSums(matrix(apply(zz, 3, rowSums), nrow = nrow(zz)))
#     cs <- cs[(1:3 - 1 + frame) %% 3 + 1]
      cs
    })
                       )
    } else if(!is.null(coordinates) & !is.null(riboDat)) {
      z <- plotCDS(coordinates, riboDat, lengths = length, plot = FALSE)
      frameCounts <- rowSums(sapply(z[[1]], function(zz) rowSums(zz[[1]])))
    } else stop("Either coordinates and riboDat parameters, or rC parameter, must be supplied to estimate the frame shift")
    frameCounts
  })

  colnames(frameCounts) <- lengths

  frameCounts <- rbind(frameCounts, frame.ML = apply(frameCounts, 2, function(x) which.max(x)) - 1)
  
  #frameShift <- which.max(frameCounts) - 1
  #c(frameShift, max(frameCounts) / sum(frameCounts))

  
  #colnames(frameShifts) <- as.character(lengths)
  #rownames(frameShifts) <- c("frameShift", "weighting")
  #frameShifts
  frameCounts
}


plotFS <- function(fS, lengths, legend.text = c("Frame 0", "Frame 1", "Frame 2"), ...)
  {
    if(!missing(lengths)) colsel <- which(colnames(fS) %in% as.character(lengths)) else colsel <- 1:ncol(fS)   
    barplot(fS[1:3, colsel], beside = TRUE, col = rainbow(3, s = 0.7), border = rainbow(3, s = 0.7), legend.text = legend.text, ...)

    
#    if(legend) legend(x = "topright", bty = "n", , fill = rainbow(3, s = 0.7), border = rainbow(3, s = 0.7), )
  }
  
