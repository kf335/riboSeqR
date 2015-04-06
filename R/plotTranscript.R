plotTranscript <-
function(transcript, coordinates, annotation, riboData, length = 27, frameShift = 0, cap, riboScale, rnaScale, xlim, main, note = "", ...) {
#  refseq <- solveMatch$refseq[solveMatch[,1] == transcript]
#  annotation <- solveMatch$annotation[solveMatch[,1] == transcript]

  alignments <- riboData@riboGR
  if(length(alignments) > 0 & is.null(dev.list()))
    par(mfrow = c(length(alignments),1), mar = c(2, 4, 2, 2))

  transceiling <- max(sapply(1:length(alignments), function(ii) {
    align <- alignments[[ii]]
    max(end(align)[as.character(seqnames(align)) == transcript])
  }))

  if(length(riboData@rnaGR) > 0) {
    transceiling <- max(c(transceiling, sapply(1:length(riboData@rnaGR), function(ii) {
      align <- riboData@rnaGR[[ii]]
      max(end(align)[as.character(seqnames(align)) == transcript])
    })))
  }
  
  if(!missing(coordinates)) transceiling <- max(c(transceiling, end(coordinates)[as.character(seqnames(coordinates)) == transcript]))
  
  for(ii in 1:length(alignments)) {    
    
    plus27GRL <- alignments[[ii]]
    allen <- plus27GRL[which(as.character(seqnames(plus27GRL)) == transcript & width(plus27GRL) == length)]
    
    createMat <- function(allen) {
      #if(any(!is.null(allen$startCodon))) startCodon <- allen$startCodon else startCodon <- 0
      #tabz <- table(start(allen) - startCodon)
      tabz <- table(start(allen) - frameShift)
      tabzz <- rep(0, ceiling(transceiling / 3) * 3)
      tabzz[as.numeric(names(tabz))] <- tabz  
      matz <- matrix(tabzz, nrow = 3);
      colnames(matz) <- (1:(ceiling(transceiling / 3) * 3))[1:(ceiling(transceiling / 3)) * 3 - 2] + 1
      matz
    }
    
    matz <- createMat(allen)
    if(!missing(cap)) matz[matz > cap] <- cap
#    if(ii < length(alignments)) colnames(matz) <- NULL
    if(missing(xlim))
      xlim <- c(0, ncol(matz) * 3)

    if(length(riboData@rnaGR) > 0) {      
      covmRNA <- riboData@rnaGR[[ii]]
      cov <- coverage(covmRNA[which(as.character(seqnames(covmRNA)) == transcript)])
      cov <- as.integer(cov[[which(names(cov) == transcript)]])      
    } else cov <- 0
    
    if(!missing(riboScale)) maxribo <- riboScale[ii] else maxribo <- max(matz)
    if(!missing(rnaScale)) maxrna <- rnaScale[ii] else maxrna <- max(cov)
    if(missing(main)) main = paste(names(alignments)[ii], " :: ", transcript, sep = "")
    ymax <- max(pretty(0:maxribo))
    
    plot(NA, NA, axes = FALSE, ylim = c(0, ymax * c(1, 1.2)[as.integer(ii == 1) + 1]), xlim = c(floor(xlim[1]), ceiling(xlim[2])), xlab = "", ylab = "")

    bp <- barplot(matz, beside = TRUE, plot = FALSE, space = c(0,1), width = 0.75, xlim = xlim)

    if(length(riboData@rnaGR) > 0) {      
      pretcov <- pretty(0:maxrna)
      cov <- cov / max(pretty(0:maxrna)) * max(pretty(0:maxribo))
      suppressWarnings(axis(side = 4, at = c(pretcov / max(pretcov)) * max(pretty(0:maxribo)), labels = pretcov, ...))
      rect(as.vector(bp)[1:length(cov)] - 0.5, 0, as.vector(bp)[1:length(cov)] + 0.5, cov, col = "light grey", border = "light grey")
    }

    bp <- barplot(matz, beside = TRUE, border = c("red", "green", "blue"), col = c("red", "green", "blue"), main = main, axes = FALSE, ylim = c(0, ymax * c(1, 1.2)[as.integer(ii == 1) + 1]), 
                  xlim = c(floor(xlim[1] / 3 * 3), ceiling(xlim[2] / 3 * 3)),
                  add = TRUE, plot = TRUE, space = c(0,1), width = 0.75)
    
    if(ii == 1 & !missing(annotation)) {
      rect(xleft = as.vector(bp)[start(annotation)], xright = as.vector(bp)[end(annotation)], ybottom = ymax *1.1, ytop = ymax * 1.15, col = "turquoise")    
      suppressWarnings(text(x = mean(as.vector(bp)[c(start(annotation), end(annotation))]), y = ymax * (1.15 + 1.1) / 2, labels = seqnames(annotation), ...))
    }
    
    suppressWarnings(axis(side = 2, pretty(0:ymax), ...))
    #if(ii == 1) text(x = 20, y = max(matz) * c(1.2, 1.1, 0.9), labels = c(refseq, annotation, note), pos = 4, cex = 0.7)
    
    if(!missing(coordinates)) {
      if(is.list(coordinates)) tw <- coordinates[[ii]][seqnames(coordinates[[ii]]) == transcript,] else tw <- coordinates[seqnames(coordinates) == transcript,]
      if(length(tw) > 0) {
        segcols <- rgb(c(0.8,0,0), c(0,0.8,0), c(0,0,0.8), alpha = 0.8)
        segments(as.vector(bp)[start(tw)] - 1, ymax * (1.015 + ((start(tw) - 1) %% 3) * 0.005), as.vector(bp)[end(tw)] - 1, ymax * (1.015 + ((start(tw) - 1) %% 3) * 0.005), col = segcols[1 + ((start(tw) - 1) %% 3)], lwd = 3)
      }
    }  
  }
  NULL
}
