frameCounting <-
function(riboDat, fastaCDS, lengths = 25:30)
  {
    message("Calling frames...", appendLF = FALSE)

    if(!("frame" %in% names(values(fastaCDS)))) fastaCDS$frame <- (start(fastaCDS) - 1) %% 3
    
    fr0GR <- fastaCDS[fastaCDS$frame == 0]
    fr1GR <- fastaCDS[fastaCDS$frame == 1]
    fr2GR <- fastaCDS[fastaCDS$frame == 2]

    frameCalls <- lapply(riboDat@riboGR, function(gral, lengths) {
#      gral$frame <- (start(gral) - 2) %% 3
      gral$frame <- (start(gral) - 1 - 0) %% 3
      if(length(fr0GR) > 0) fr0 <- .tableOverlaps(fr0GR, gral, lengthRange = lengths) else fr0 <- new("riboCoding", hits = array(dim=c(0,1,3,length(lengths))), unqHits = array(dim=c(0,1,3,length(lengths))))
      message(".", appendLF = FALSE)
      gral$frame <- (start(gral) - 1 - 1) %% 3
      if(length(fr1GR) > 0) fr1 <- .tableOverlaps(fr1GR, gral, lengthRange = lengths) else fr1 <- new("riboCoding", hits = array(dim=c(0,1,3,length(lengths))), unqHits = array(dim=c(0,1,3,length(lengths))))
      message(".", appendLF = FALSE)
      gral$frame <- (start(gral) - 1 - 2) %% 3
      if(length(fr2GR) > 0) fr2 <- .tableOverlaps(fr2GR, gral, lengthRange = lengths) else fr2 <- new("riboCoding", hits = array(dim=c(0,1,3,length(lengths))), unqHits = array(dim=c(0,1,3,length(lengths))))
      message(".", appendLF = FALSE)
      rc <- new("riboCoding",
          CDS = c(fr0@CDS, fr1@CDS, fr2@CDS),
          hits = do.call("abind", args = list(list(fr0@hits, fr1@hits, fr2@hits), along = 1)),
          unqHits = do.call("abind", args = list(list(fr0@unqHits, fr1@unqHits, fr2@unqHits), along = 1)))
      rc
    }, lengths = lengths)
    message(".done!", appendLF = TRUE)
    fCs <- new("riboCoding")
    fCs@hits <- do.call("abind", args = list(lapply(frameCalls, function(x) x@hits), along = 2))
    fCs@unqHits <- do.call("abind", args = list(lapply(frameCalls, function(x) x@unqHits), along = 2))
    fCs@CDS <- frameCalls[[1]]@CDS
    fCs@replicates <- riboDat@replicates
    fCs
  }
