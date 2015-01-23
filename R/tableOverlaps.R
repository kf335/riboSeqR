.tableOverlaps <-
function(grf, gral, uniqueHits = FALSE, lengthRange = 25:30) {
  #gral.f <- gral[gral$frame == frame]
  gral <- gral[width(gral) %in% lengthRange]
  spl27.f <- split(start(gral), factor(as.character(seqnames(gral)), levels = as.character(seqnames(grf))))
  spl27.frame <- split(gral$frame, factor(as.character(seqnames(gral)), levels = as.character(seqnames(grf))))
  spl27.length <- split(width(gral), factor(as.character(seqnames(gral)), levels = as.character(seqnames(grf))))
  
  splfr0 <- split(start(grf), factor(as.character(seqnames(grf)), levels = as.character(seqnames(grf))))
  splfr0e <- split(end(grf), factor(as.character(seqnames(grf)), levels = as.character(seqnames(grf))))  
  tov <- lapply(1:length(spl27.f), function(ii) {
    if(length(splfr0[[ii]]) > 0) {      
      z <- data.frame(startCDS = findInterval(spl27.f[[ii]], splfr0[[ii]]), endCDS = findInterval(spl27.f[[ii]], splfr0e[[ii]]), start = spl27.f[[ii]], frame = spl27.frame[[ii]], length = spl27.length[[ii]])
      z <- z[z[,1] != 0 & z[,1] == z[,2] + 1,]
      lfTab <- as.data.frame(
                 table(factor(z[,1], levels = 1:length(splfr0[[ii]])), factor(z$frame, levels = 0:2), factor(z$length, levels = lengthRange))
                 )
      lfValues <- matrix(lfTab[,4], nrow = length(splfr0[[ii]]))
      colnames(lfValues) <- paste(lfTab[,2], lfTab[,3], sep = ":")[lfTab[,1] == 1]      
      
      zunq <- z[!duplicated(z[,c("start", "length")]),]
      lfTabUnq <- as.data.frame(table(factor(zunq[,1], levels = 1:length(splfr0[[ii]])), factor(zunq$frame, levels = 0:2), factor(zunq$length, levels = lengthRange)))
      lfValuesUnq <- matrix(lfTabUnq[,4], nrow = length(splfr0[[ii]]))
      colnames(lfValuesUnq) <- paste(lfTabUnq[,2], lfTabUnq[,3], sep = ":")[lfTabUnq[,1] == 1]      
      
#      }
      list(lfValues = lfValues, lfValuesUnq = lfValuesUnq)
    } else NULL
  })

  hits <- do.call("rbind", (lapply(tov, function(x) x$lfValues)))
  hits <- array(hits, dim = c(nrow(hits), 1, 3, length(lengthRange)), dimnames = list(NULL, NULL, 0:2, lengthRange))
  
  unqHits <- do.call("rbind", (lapply(tov, function(x) x$lfValuesUnq)))
  unqHits <- array(hits, dim = c(nrow(hits), 1, 3, length(lengthRange)), dimnames = list(NULL, NULL, 0:2, lengthRange))
  message(".", appendLF = FALSE)  

  new("riboCoding", CDS = grf, hits = hits, unqHits = unqHits)
}
