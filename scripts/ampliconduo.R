library(AmpliconDuo)
library(ggplot2)
library(xtable)
library(data.table)

# apparently ggplot opens a device to get for example hight and width for plots if not specified
# which in turn does produce an additional Rplots.pdf in the working dir. 
# The following line should suppress this.
pdf(NULL)

log <- file(toString(snakemake@log), open="wt")
sink(log, append = TRUE)
sink(log, type="message", append =TRUE)

# Script for statistical analysis of the filtering process (filtered vs unfiltered),
# it will also use Fisher's exact test to find significantly deviating read numbers
# between split-samples, which could be used to implement an additional filtering step.
unfiltered.table <- fread(snakemake@input[["unfiltered_table"]],
                          stringsAsFactors=FALSE, data.table=FALSE)
filtered.table <- fread(snakemake@input[["filtered_table"]],
                        stringsAsFactors=FALSE, data.table=FALSE)
figure.folder <- unlist(strsplit(toString(snakemake@output),
                                 split="/AmpliconDuo.RData"))

amp.duo <- function(table, figure.folder, saving.format, file.name, p.corr) {
  intable <- table[, -1]
  names <- names(intable)
  names <- names[seq(1, length(names), by=2)]
  names <- sapply(X=names, FUN=function(x){
  unlist(strsplit(x, split="+.A$"))
  })
  splitsamples <-  lapply(names, function(x) grep(x, colnames(intable), value=T))
  allzero <- unlist(lapply(splitsamples, function(x) all(rowSums(intable[,x])==0)))
  if(any(allzero)){
    print("These splitsamples do not have any OTUS/ASVs in common:")
    print(splitsamples[allzero])
    intable <- intable[,unlist(splitsamples[!allzero])]
    names <- names[!allzero]
  }
  ampliconduos <- ampliconduo(intable, sample.names=names,
                              correction=p.corr)
  if (snakemake@params[["plot_ampduo"]]) {
    nrow = as.integer(sqrt(length(names)))
    plotAmpliconduo.set(ampliconduos, nrow=nrow, save=TRUE,
                        path=figure.folder, format=saving.format,
                        file.name=file.name)
  }
  return(ampliconduos)
}

discordance <- function(ampliconduos, figure.folder, file.name) {
  discordance.delta(ampliconduos, corrected=TRUE, printToTex=TRUE,
                    directory=figure.folder, file.name=file.name)
}

amp.duo.unfiltered = amp.duo(unfiltered.table, figure.folder,
                             snakemake@params[["saving_format"]],
                             paste("ampliconduo_unfiltered_",
                                   Sys.Date(), sep=""),
                             snakemake@params[["p_corr"]])

discordance(amp.duo.unfiltered, figure.folder,
            paste("discordanceDelta_unfiltered_", Sys.Date(), sep=""))

amp.duo.filtered = amp.duo(filtered.table, figure.folder,
                           snakemake@params[["saving_format"]],
                           paste("ampliconduo_filtered_", Sys.Date(), sep=""),
                           snakemake@params[["p_corr"]])

discordance(amp.duo.filtered, figure.folder,
            paste("discordanceDelta_filtered_", Sys.Date(), sep=""))

save(amp.duo.unfiltered, amp.duo.filtered, file=paste(figure.folder,
                                                        "/AmpliconDuo.RData",
                                                        sep=""))

sink()
sink(type="message")
