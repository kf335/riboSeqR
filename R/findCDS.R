findCDS <-
function(fastaFile, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))
  {
    fasta <- scan(fastaFile, sep = "\n", what = "character")
    fastaHeads <- gsub(">", "", fasta[grep(">", fasta)])
    headloc <- grep(">", fasta)
    fastaSeq <- sapply(1:length(headloc), function(ii) toupper(paste(fasta[(headloc[ii] + 1):c(headloc[-1] - 1, length(fasta))[ii]], collapse = "")))
    
    codpos <- paste("(", paste(startCodon, collapse = "|"), " )(.{3} )*?(", paste(stopCodon, collapse = "|"), "|.{0,2}$)", sep = "")
    frame0 <- gregexpr(codpos, gsub("(.{3})", "\\1 ", fastaSeq))
    frame1 <- gregexpr(codpos, gsub("(.{3})", "\\1 ", substr(fastaSeq, 2, nchar(fastaSeq))))
    frame2 <- gregexpr(codpos, gsub("(.{3})", "\\1 ", substr(fastaSeq, 3, nchar(fastaSeq))))
    
    fr0mat <- lapply(frame0, function(x) if(x[1] != -1) cbind((as.integer(x) + 3) / 4 * 3 -2, (attr(x, "match.length") + 1) / 4 * 3) else matrix(NA, ncol = 2, nrow = 0))
    fr0GR <- GRanges(seqnames = factor(rep(fastaHeads, sapply(fr0mat, nrow)), levels = fastaHeads),
                     IRanges(start = do.call("rbind", fr0mat)[,1], width = do.call("rbind", fr0mat)[,2]),
                     frame = rep(0, sum(sapply(fr0mat, nrow))))
    fr1mat <- lapply(frame1, function(x) if(x[1] != -1) cbind((as.integer(x) + 3) / 4 * 3 - 1, (attr(x, "match.length") + 1) / 4 * 3) else matrix(NA, ncol = 2, nrow = 0))
    fr1GR <- GRanges(seqnames = factor(rep(fastaHeads, sapply(fr1mat, nrow)), levels = fastaHeads),
                     IRanges(start = do.call("rbind", fr1mat)[,1], width = do.call("rbind", fr1mat)[,2]),
                     frame = rep(1, sum(sapply(fr1mat, nrow))))
    fr2mat <- lapply(frame2, function(x) if(x[1] != -1) cbind((as.integer(x) + 3) / 4 * 3, (attr(x, "match.length") + 1) / 4 * 3) else matrix(NA, ncol = 2, nrow = 0))
    fr2GR <- GRanges(seqnames = factor(rep(fastaHeads, sapply(fr2mat, nrow)), levels = fastaHeads),
                     IRanges(start = do.call("rbind", fr2mat)[,1], width = do.call("rbind", fr2mat)[,2]),
                     frame = rep(2, sum(sapply(fr2mat, nrow))))

    gr <- c(fr0GR, fr1GR, fr2GR)
    gr <- gr[order(as.integer(seqnames(gr)), start(gr), end(gr))]
    gr
  }
