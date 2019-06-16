library("data.table")
library("tidyverse")

ncbi_get_tax <- function() {
    ncbi_tax_dir <- dirname(snakemake@params[["db_path"]])
    blast_table <- fread(snakemake@input[[1]])
    ranked_lineage_file_index  <- grep("rankedlineage.dmp" , list.files(ncbi_tax_dir , full.names = T))

    ranked_lineage_file <- list.files(ncbi_tax_dir , full.names = T)[ranked_lineage_file_index]

    ranked_lineage <- fread(ranked_lineage_file , data.table = F, na.strings = "") %>%
      as_tibble() %>%
      select(seq(from = 1, to = 20 , by = 2)) %>%
      `colnames<-`(c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"))

    ncbi_tax_ids <- as.numeric(ncbi_tax_ids) %>% as_tibble() %>% `colnames<-`(c("query_tax_id"))
    mapped <- ranked_lineage %>%
      right_join(blast_table , by = c("tax_id" = "taxonomy"))

    d <- unite(mapped, "taxonomy", rev(c(colnames(mapped[2:10]))), remove=TRUE, sep = ";")
    df <- subset(d, select=c(3:14,2,15:length(d)))
    df$taxonomy <- gsub("NA;", "", df$taxonomy)

    fwrite(df, file = snakemake@output[[1]], na="")
}

if (snakemake@params[['db']] == 'NCBI') {
    ncbi_get_tax()
} else {
    file.copy(snakemake@input[[1]], snakemake@output[[1]])
}
