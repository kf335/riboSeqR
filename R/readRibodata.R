.readAlignments <-
  function(filenames, columns, header, seqnames, zeroIndexed) {
    positiveOnly = TRUE
    if(missing(seqnames)) seqnames <- NULL
    GRLs <- lapply(filenames, function(file) {
      fdat <- read.delim(file, as.is = TRUE, header = header)[,columns]
      if(positiveOnly) fdat <- fdat[fdat[,1] == "+",]
      if(!is.null(seqnames)) {
        fdat <- fdat[fdat[,2] %in% seqnames,]
        fsn <- seqnames
      } else fsn <- unique(fdat[,2])
      
      GRL <- GRanges(                 
               seqnames = factor(fdat[,2], levels = fsn),
               IRanges(start = fdat[,3] + as.integer(zeroIndexed), width = nchar(fdat[,4])),
               strand = fdat[,1])
      
      message(".", appendLF = FALSE)
      GRL
    })
    names(GRLs) <- gsub("^.*/", "", filenames)
    GRLs
  }

readRibodata <-
function(riboFiles, rnaFiles, columns = c(strand = 1, seqname = 2, start = 3, sequence = 4), zeroIndexed = TRUE, header = FALSE, replicates, seqnames) {
  riboDat <- new("riboData", replicates = as.factor(replicates))
  message("Reading ribosomal files...", appendLF = FALSE)
  riboDat@riboGR <- .readAlignments(riboFiles, columns = columns, header = header, seqnames, zeroIndexed = zeroIndexed)
  message("done!")
  if(!missing(rnaFiles)) {
    message("Reading rna files...", appendLF = FALSE)
    riboDat@rnaGR <- .readAlignments(rnaFiles, columns = columns, header = header, seqnames, zeroIndexed = zeroIndexed)
    message("done!")
  }
  riboDat
}
