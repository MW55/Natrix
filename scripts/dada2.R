library(dada2); packageVersion("dada2")
library(ShortRead)

sink(toString(snakemake@log), append = TRUE)
sample.names <- sapply(strsplit(basename(snakemake@input[["forward"]]), "_[12]_cut.fastq"), `[`, 1)
fnFs <- snakemake@input[["forward"]]

names(fnFs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(fnFs, nbases=1e8, multithread=TRUE, randomize = TRUE, verbose = TRUE)

if(snakemake@params[["paired_end"]] == TRUE){
  fnRs <- snakemake@input[["reverse"]]
  # Learn reverse error rates 
  names(fnRs) <- sample.names
  errR <- learnErrors(fnRs, nbases=1e8, multithread=TRUE, randomize = TRUE, verbose = TRUE)
}

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(fnFs[[sam]], verbose = TRUE)
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose = TRUE)
  if (snakemake@params[["paired_end"]] == TRUE) {
    derepR <- derepFastq(fnRs[[sam]], verbose = TRUE)
    ddR <- dada(derepR, err=errR, multithread=TRUE, verbose = TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
  } else {
    mergers[[sam]] <- ddF
  }
}
rm(derepF)
if (snakemake@params[["paired_end"]] == TRUE) {
  rm(derepR)
}
seqtab <- makeSequenceTable(mergers)

for (i in seq(nrow(seqtab))) {
  sample_df <- t(as.data.frame(seqtab[i,]))
  only_present_seqs <- t(as.data.frame(sample_df[,-(which(colSums(sample_df)==0))]))
  seq_uniq <- getUniques(only_present_seqs)
  seq_names <- paste0(as.character(seq(length(only_present_seqs))), ";size=", seq_uniq, ";")
  uniquesToFasta(seq_uniq, fout = snakemake@output[[i]], ids = seq_names)
}
sink()