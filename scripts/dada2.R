library(dada2); packageVersion("dada2")
library(ShortRead)

sink(snakemake@log, append = TRUE)
sample.names <- sapply(strsplit(basename(snakemake@input[["forward"]]), "_[12]_cut.fastq"), `[`, 1)

#plotQualityProfile(fnFs[1:2])
derepF1 <- derepFastq(snakemake@input[["forward"]], verbose=TRUE)
errF <- learnErrors(derepF1, multithread=TRUE, verbose=TRUE)
dadaFs <- dada(derepF1, err=errF, multithread=TRUE, verbose=TRUE)
print("DADA2 forward reads:")
dadaFs[[1]]

if(snakemake@params[["paired_end"]] == TRUE){
  derepF2 <- derepFastq(snakemake@input[["reverse"]], verbose=TRUE)
  errR <- learnErrors(derepF2, multithread=TRUE, verbose=TRUE)
  dadaRs <- dada(derepF2, err=errR, multithread=TRUE, verbose=TRUE)
  mergers <- mergePairs(dadaFs, snakemake@input[["forward"]], dadaRs, snakemake@input[["reverse"]], verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
  print("DADA2 reverse reads:\n")
  dadaRs[[1]]
} else {
  seqtab <- makeSequenceTable(dadaFs)
}

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

for (i in seq(nrow(seqtab))) {
  sample_df <- t(as.data.frame(seqtab[i,]))
  only_present_seqs <- t(as.data.frame(sample_df[,-(which(colSums(sample_df)==0))]))
  seq_uniq <- getUniques(only_present_seqs)
  seq_names <- paste0(as.character(seq(length(only_present_seqs))), ";size=", seq_uniq, ";")
  uniquesToFasta(seq_uniq, fout = snakemake@output[[i]], ids = seq_names)
  
}
  