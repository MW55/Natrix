library(dada2); packageVersion("dada2")
library(ShortRead)

# WIP

#pathf<- "~/projects/Natrix/results/assembly/*/*_[AB]_[1]_cut.fastq"
#pathr <- "~/projects/Natrix/results/assembly/*/*_[AB]_[2]_cut.fastq"

#fnFs <- sort(Sys.glob(pathf))
#fnRs <- sort(Sys.glob(pathr))

print(snakemake@input[["forward"]])
print(snakemake@output)

sample.names <- sapply(strsplit(basename(snakemake@input[["forward"]]), "_[12]_cut.fastq"), `[`, 1)

#plotQualityProfile(fnFs[1:2])

derepF1 <- derepFastq(snakemake@input[["forward"]], verbose=TRUE)
derepF2 <- derepFastq(snakemake@input[["reverse"]], verbose=TRUE)

errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepF2, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(derepF1, err=errF, multithread=TRUE)
dadaRs <- dada(derepF2, err=errR, multithread=TRUE)

# Few successful merges, might be better to use PandaSeq 
mergers <- mergePairs(dadaFs, snakemake@input[["forward"]], dadaRs, snakemake@input[["reverse"]], verbose=TRUE)

#seqtabF <- t(makeSequenceTable(dadaFs))
#seqtabR <- t(makeSequenceTable(dadaRs))

seqtab <- makeSequenceTable(mergers)

for (i in seq(nrow(seqtab))) {
  sample_df <- t(as.data.frame(seqtab[i,]))
  only_present_seqs <- t(as.data.frame(sample_df[,-(which(colSums(sample_df)==0))]))
  seq_uniq <- getUniques(only_present_seqs)
  seq_names <- paste0(as.character(seq(length(only_present_seqs))), ";size=", seq_uniq, ";")
  uniquesToFasta(seq_uniq, fout = snakemake@output[i], ids = seq_names)
  
}
  