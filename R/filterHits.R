filterHits <-
  function(fCs, lengths = 27, frames, hitMean = 10, unqhitMean = 1, ratioCheck = TRUE, fS) {
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
    ffCs <- fCs[filtration,]
    
    if(ratioCheck) {
      trues <- sapply(1:length(lengths), function(ll, fS) {
        len <- lengths[ll]
        lenfram <- frames[[ll]]
        keep <- rep(TRUE, nrow(ffCs@hits))
        if(length(lenfram) != 3 || !all(sort(lenfram) == 0:2)) {
          sumhits <- apply(ffCs@hits[,,,as.character(len),drop = FALSE], c(1,3), sum)
          fps <- which(rowSums(sumhits[,as.character(lenfram),drop = FALSE]) < apply(sumhits[,!colnames(sumhits) %in% as.character(lenfram),drop = FALSE], 1, max))
          if(length(fps) > 0) {
            reject <- fps[sapply(fps, function(ff, fS) {
              maxframe <- which.max(sumhits[ff,]) - 1
              maxfS <- which.max(fS[,as.character(len)]) - 1
              rmap <- sumhits[ff, 1 + ((0:2 - (maxfS - maxframe)) %% 3)]
              compdat <- rbind(c(rmap[maxfS + 1], sum(sumhits[ff,as.character(lenfram)])),
                               c(max(fS[,as.character(len)]), sum(fS[1 + ((lenfram - (maxframe - maxfS)) %% 3), as.character(len)])))
              if(which.max(compdat[,2] / compdat[,1]) == 2 | chisq.test(compdat)$p.value > 0.05) return(TRUE) else return(FALSE)
            }, fS = fS)]
            if(length(reject) > 0) keep[reject] <- FALSE
          }
        }
        return(keep)
      }, fS = fS)
      ffCs <- ffCs[rowSums(trues) > 0,]
    }
    
    ffCs
  }
  
