library(dada2)
library(ShortRead)

log <- file(toString(snakemake@log), open="wt")
sink(log, append = TRUE)
sink(log, type="message", append =TRUE)

paired.end <- snakemake@params[["paired_end"]]
min.overlap <- snakemake@params[["minoverlap"]]
split.samples <- snakemake@params[["splitsamples"]]=="split_sample"
sample.names <- sapply(strsplit(basename(snakemake@input[["fwd"]]), "_[12]_cut.fastq"), `[`, 1)
fnFs <- snakemake@input[["fwd"]]
names(fnFs) <- sample.names

set.seed(100)
# Filtering was done before with prinseq and cutadapt
# Learn forward error rates
errF <- learnErrors(fnFs, nbases=1e8, multithread=TRUE, randomize = TRUE, verbose = TRUE)

if(paired.end){
  fnRs <- snakemake@input[["rev"]]
  # Learn reverse error rates
  names(fnRs) <- sample.names
  errR <- learnErrors(fnRs, nbases=1e8, multithread=TRUE, randomize = TRUE, verbose = TRUE)
}

if(split.samples){
  samples <- unlist(unique(strsplit(sample.names, "_[AB]$")))
  sample.names2 <- lapply(samples, function(x) grep(paste("^",x,sep=""), sample.names, value=TRUE))
  names(sample.names2) <- samples
} else {
  sample.names2 <- sample.names
  names(sample.names2) <- sample.names
}

n <- 1
for(i in 1:length(sample.names2)) {
  result <- c()
  sam <- sample.names2[i]
  cat("Processing:", names(sam), "\n")
  # dereplication done on the fly since dada 1.12
  dadaF <- dada(fnFs[unlist(sam)], err=errF, multithread=TRUE, verbose = TRUE)
  print("Forward After DADA2 dereplication and denoising:")
  print(dadaF)
  if (paired.end == TRUE) {
    dadaR <- dada(fnRs[unlist(sam)], err=errR, multithread=TRUE, verbose = TRUE)
    print("Reverse After DADA2 dereplication and denoising:")
    print(dadaR)
    result <- mergePairs(dadaF, fnFs[unlist(sam)], dadaR, fnRs[unlist(sam)], verbose = TRUE, minOverlap=min.overlap)
    if(!is(result, "list")){result <- list(result)}
  } else {
    if(is(dadaF, "dada")){dadaF <- list(dadaF)}
    seqtab <- lapply(dadaF, function(x) t(makeSequenceTable(x)))
    result <- lapply(seqtab, function(x) data.frame(abundance = x[,1], sequence=rownames(x)))
  }
  print("Sequences left after assembly:")
  lapply(result, function(x) print(nrow(x)))
  for(j in 1:length(result)){
    if(nrow(result[[j]]) == 0)
    {
      print(paste0("No sequences left for sample ", names(sam)))
      cat(NULL, file=snakemake@output[[j]])
    } else {
      seq_names = paste0(seq(nrow(result[[j]])), ";size=", result[[j]]$abundance, ";")
      uniquesToFasta(result[[j]], fout = snakemake@output[[n]], ids = seq_names)
    }
    n <- n+1
  }
}

sink()
sink(type="message")
