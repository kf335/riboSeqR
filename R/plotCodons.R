plotCDS <- function(coordinates, riboDat, lengths = 27, min5p = -20, max5p = 200, min3p = -200, max3p = 20, cap, main = "", plot = TRUE, ...) {

  min5p <- floor(min5p / 3) * 3; max5p <- (ceiling(max5p / 3) * 3) - 1
  min3p <- floor(min3p / 3) * 3; max3p <- (ceiling(max3p / 3) * 3) - 1

  alignments <- riboDat@riboGR
  if(plot && length(alignments) > 0 & is.null(dev.list()))
    par(mfrow = c(length(alignments), 1))
  if(missing(cap)) cap <- NULL
  
  lenmatz <- lapply(lengths, function(length) {
    alignments <- lapply(alignments, function(al) al[width(al) == length])  
    matzes <- lapply(1:length(alignments), function(ii) {
      allen <- alignments[[ii]]
      redcod <- coordinates[as.character(seqnames(coordinates)) %in% seqlevels(allen),]
#      allen <- allen[match(seqnames(redcod), seqlevels(allen))]
      starts <- split(start(allen), as.character(seqnames(allen)))
      rcmat <- match(names(starts), as.character(seqnames(redcod)))

      strod <- vector("list", length(redcod))
      strod[rcmat[!is.na(rcmat)]] <- starts[!is.na(rcmat)]
      
      
      makeMatz <- function(whichEnd, min, max) {
        if(whichEnd == "5p") fixedPoints <- start(redcod) else fixedPoints <- end(redcod) - 2
        widthcod <- width(redcod)
#        lencods <- cbind(unlist(strod) - rep(fixedPoints, sapply(strod, length)),
#                         1 / (rep(sapply(strod, length),
#                                  sapply(strod, length)) / rep(width(redcod), sapply(strod, length))), rep(1:length(strod), sapply(strod, length)))

        adjstrod <- split(unlist(strod) - rep(fixedPoints, sapply(strod, length)), factor(rep(1:length(strod), sapply(strod, length)), levels = 1:length(strod)))
        ctstrod <- do.call("rbind", lapply(adjstrod, function(x) if(is.null(x)) return(rep(0, max - min + 1)) else table(factor(x[x >= min & x<= max], levels = min:max))))
        weights <- 1 / (sapply(strod, length) / (rep(width(redcod)/1000)))
        weightMean <- apply(ctstrod[weights < Inf,,drop = FALSE], 2, weighted.mean, weights[weights < Inf]) 
                                                                          
#        lencods <- lencods[lencods[,1] >= min & lencods[,1] <= max,]
#        spllen <- split(lencods[,2], factor(lencods[,1], levels = min:max))
#        weightsum <- sapply(spllen, function(x) sum(1/x))    
        matz <- matrix(weightMean, nrow = 3);
        colnames(matz) <- (min:max)[1:((max - min + 1) / 3) * 3 - 2]
        matz
      }
      
      if(min5p < max5p) matz5p <- makeMatz("5p", min = min5p, max = max5p) else matz5p <- NULL
      if(min3p < max3p) matz3p <- makeMatz("3p", min = min3p, max = max3p) else matz3p <- NULL
      matz <- cbind(matz5p, matrix(NA, ncol = max(1, min(25, round(sum(ncol(matz5p), ncol(matz5p)) * 0.1))), nrow = 3), matz3p)      
      if(!is.null(cap)) matz[matz > cap] <- cap
      
      if(plot) {
        bw <- barplot(matz, beside = TRUE, col = rainbow(3, s = 0.7), border = rainbow(3, s = 0.7), main = main, ylim = c(0, max(matz, na.rm = TRUE)), axisnames = FALSE, space = c(0, 0.25), ylab = "Mean number of reads", xlab = "Base position relative to CDS",...)        
        if(min5p < max5p) {
          pretty5p <- pretty(min5p:max5p); pretty5p <- pretty5p[pretty5p >= min5p]; pretty5p <- pretty5p[pretty5p <= max5p]; pretty5p <- unique(pretty5p)
          prettystart <- pretty5p; prettystart[prettystart == 0] <- "start"
          axis(side = 1, at = as.vector(bw)[(pretty5p - min5p + 1)], labels = prettystart, lwd = 0, lwd.ticks = 1, ...)
        }
        if(min3p < max3p) {
          pretty3p <- pretty(min3p:max3p); pretty3p <- pretty3p[pretty3p >= min3p]; pretty3p <- pretty3p[pretty3p <= max3p]; pretty3p <- unique(pretty3p)
          prettystop <- pretty3p; prettystop[prettystop == 0] <- "stop"
          axis(side = 1, at = as.vector(bw[,max(which(is.na(matz[1,])) + 1):ncol(bw)])[pretty3p - min3p + 1], labels = prettystop, lwd = 0, lwd.ticks = 1, ...)
        }
        
      }
      matz
    })
  })
  invisible(lenmatz)
}
