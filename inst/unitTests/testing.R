
test_findCDS <- function() {
  fcds <- findCDS("../extdata/unit_test.fasta", startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))
  checkTrue(all(as.data.frame(fcds) ==
      data.frame(seqnames = factor(c("ut_1", "ut_3"), levels = c("ut_1", "ut_2", "ut_3")), start = c(1,2), end = c(7,10), width = c(7,9), strand = factor(c("*", "*"), levels = c("+", "-", "*")), frame = c(0, 1))))
}

test_readRibodata <- function() {
  datadir <- system.file("extdata", package = "riboSeq")
  chlamyFasta <- paste(datadir, "/rsem_chlamy236_deNovo.transcripts.fa", sep = "")  
  fastaCDS <- findCDS(fastaFile = chlamyFasta,
                      startCodon = c("ATG"),
                      stopCodon = c("TAG", "TAA", "TGA"))
  rR <- readRibodata(paste(datadir, "/chlamy236_plus_deNovo_plusOnly_Index17_short", sep = ""), replicates = c(1))
  fCs <- frameCounting(rR, fastaCDS)
  checkTrue(all(frameShift(rC = fCs)[,2] == c(0,2,0,1)))
}
