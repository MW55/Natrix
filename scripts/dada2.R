library(dada2)

# WIP

pathf<- "~/projects/Natrix/results/assembly/*/*_[AB]_[1]_cut.fastq"
pathr <- "~/projects/Natrix/results/assembly/*/*_[AB]_[2]_cut.fastq"

fnFs <- sort(Sys.glob(pathf))
fnRs <- sort(Sys.glob(pathr))

sample.names <- sapply(strsplit(basename(fnFs), "_[12]_cut.fastq"), `[`, 1)

plotQualityProfile(fnFs[1:2])

derepF1 <- derepFastq(fnFs, verbose=TRUE)
derepF2 <- derepFastq(fnRs, verbose=TRUE)

errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepF2, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(derepF1, err=errF, multithread=TRUE)
dadaRs <- dada(derepF2, err=errR, multithread=TRUE)

# Few successful merges, might be better to use PandaSeq 
mergers <- mergePairs(dadaFs, fnFs, dadaRs, fnRs, verbose=TRUE)

#seqtabF <- t(makeSequenceTable(dadaFs))
#seqtabR <- t(makeSequenceTable(dadaRs))

seqtab <- t(makeSequenceTable(mergers))
