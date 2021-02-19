library(dada2)
library(ShortRead)

log <- file(toString(snakemake@log), open="wt")
sink(log, append = TRUE)
sink(log, type="message", append =TRUE)

paired_end <- snakemake@params[["paired_end"]]
min.overlap <- snakemake@params[["minoverlap"]]
sample.names <- sapply(strsplit(basename(snakemake@input[["fwd"]]), "_[12]_cut.fastq"), `[`, 1)
fnFs <- snakemake@input[["fwd"]]
names(fnFs) <- sample.names

set.seed(100)
# Filtering was done before with prinseq and cutadapt
# Learn forward error rates
errF <- learnErrors(fnFs, nbases=1e8, multithread=TRUE, randomize = TRUE, verbose = TRUE)

if(paired_end == TRUE){
  fnRs <- snakemake@input[["rev"]]
  # Learn reverse error rates
  names(fnRs) <- sample.names
  errR <- learnErrors(fnRs, nbases=1e8, multithread=TRUE, randomize = TRUE, verbose = TRUE)
}

for(i in 1:length(sample.names)) {
  result <- c()
  sam <- sample.names[i]
  cat("Processing:", sam, "\n")
  # dereplication done on the fly since dada 1.12
  dadaF <- dada(fnFs[sam], err=errF, multithread=TRUE, verbose = TRUE)
  print("Forward After DADA2 dereplication and denoising:")
  print(dadaF)
  if (paired_end == TRUE) {
    dadaR <- dada(fnRs[sam], err=errR, multithread=TRUE, verbose = TRUE)
    print("Reverse After DADA2 dereplication and denoising:")
    print(dadaR)
    result <- mergePairs(dadaF, fnFs[sam], dadaR, fnRs[sam], verbose = TRUE, minOverlap=min.overlap)
  } else {
    seqtab <- t(makeSequenceTable(dadaF))
    result <- data.frame(abundance = seqtab[,1], sequence=rownames(seqtab))
  }
  print("Sequences left after assembly:")
  print(nrow(result))
  if(nrow(result) == 0)
  {
    print(paste0("No sequences left for sample ", sam))
    cat(NULL, file=snakemake@output[[i]])
  } else {
    seq_names = paste0(seq(nrow(result)), ";size=", result$abundance, ";")
    uniquesToFasta(result, fout = snakemake@output[[i]], ids = seq_names)
  }
}

sink()
sink(type="message")
