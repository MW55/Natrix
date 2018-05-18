library(ggplot2)
library(xtable)
library(data.table)

source('bin/AmpliconDuo/R/ampliconduo.R')
source('bin/AmpliconDuo/R/discordance.delta.R')
source('bin/AmpliconDuo/R/plotAmpliconduo.set.R')

unfiltered.table <- fread(snakemake@input[['unfiltered_table']],
                          stringsAsFactors=FALSE, data.table=FALSE)
filtered.table <- fread(snakemake@input[['filtered_table']],
                        stringsAsFactors=FALSE, data.table=FALSE)
figure.folder <- unlist(strsplit(toString(snakemake@output),
                                 split = '/AmpliconDuo.RData'))

amp.duo <- function(table, figure.folder, saving.format, file.name, p.corr) {
  intable <- table[, -1]
  names <- names(intable)
  names <- names[seq(1, length(names), by = 2)]
  names <- sapply(X = names, FUN = function(x){
  unlist(strsplit(x, split = "+.A$"))
  })
  ampliconduos <- ampliconduo(intable, sample.names = names,
                              correction = p.corr)
  if (snakemake@params[['plot_ampduo']]) {
    nrow = as.integer(sqrt(length(names)))
    plotAmpliconduo.set(ampliconduos, nrow = nrow, save = TRUE,
                        path = figure.folder, format = saving.format,
                        file.name = file.name)
  }
  return(ampliconduos)
}

## discordance Delta befor filtering
discordance <- function(ampliconduos, figure.folder, file.name) {
  discordance.delta(ampliconduos, corrected = TRUE, printToTex= TRUE,
                    directory = figure.folder, file.name = file.name)
}

amp.duo.unfiltered = amp.duo(unfiltered.table, figure.folder,
                             snakemake@params[['saving_format']],
                             paste("ampliconduo_unfiltered_",
                                   Sys.Date(), sep = ""),
                             snakemake@params[['p_corr']])

discordance(amp.duo.unfiltered, figure.folder,
            paste("discordanceDelta_unfiltered_", Sys.Date(), sep = ""))

amp.duo.filtered = amp.duo(filtered.table, figure.folder,
                           snakemake@params[['saving_format']],
                           paste("ampliconduo_filtered_", Sys.Date(), sep = ""),
                           snakemake@params[['p_corr']])

discordance(amp.duo.filtered, figure.folder,
            paste("discordanceDelta_filtered_", Sys.Date(), sep = ""))

save(amp.duo.unfiltered, amp.duo.filtered, file = paste(figure.folder,
                                                        '/AmpliconDuo.RData',
                                                        sep = ''))
