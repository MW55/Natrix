# Script to indentify the primer sequences depending on the length of th primers (offset) or the sequence of the primers.
primer  <- read.table(snakemake@input[['primer_table']], header=T, stringsAsFactors=F, sep= ',')
comb.f = c() ## forward primer
comb.r = c() ## reverse primer
## pandaseq just takes the length of the primers but not the sequence, in case primer are different
## and can' be expressed by common 1L code
if(snakemake@config[['qc']][['primer_offset']] == TRUE){
  if (snakemake@config[['qc']][['mid_rm']]){ ## only forward primer
    for(i in 1:length(primer[,1])){
      tmp = primer[i,"specific_forward_primer"]
      tmp = unlist(strsplit(tmp, split = ""))
      comb.f <- c(comb.f, length(tmp))
      tmp = paste0(primer[i, "poly_N_rev"], primer[i, "specific_reverse_primer"], collapse = "")
      tmp = unlist(strsplit(tmp, split = ""))
      comb.r <- c(comb.r, length(tmp))
    }
  }else{ ### takes primer sequence
    for(i in 1:length(primer[,1])){
      tmp = paste0(primer[i,"poly_N"], primer[i,"MID"],primer[i,"specific_forward_primer"], collapse = "")
      tmp = unlist(strsplit(tmp, split = ""))
      comb.f <- c(comb.f, length(tmp))
      tmp = paste0(primer[i, "poly_N_rev"], primer[i, "specific_reverse_primer"], collapse = "")
      tmp = unlist(strsplit(tmp, split = ""))
      comb.r <- c(comb.r, length(tmp))
    }
  }
  primer <- cbind(primer, f_primer = comb.f, r_primer = comb.r)
}else{ # work with sequences not offset
  if (snakemake@config[['qc']][['mid_rm']]){
    for(i in 1:length(primer[,1])){
      comb.f <- c(comb.f,  as.character(primer[i,"specific_forward_primer"]))
      comb.r <- c(comb.r, as.character(primer[i, "specific_reverse_primer"]))
    }
  }else{
    for(i in 1:length(primer[,1])){
      comb.f <- c(comb.f, paste( primer[i,"MID"], primer[i,"specific_forward_primer"], sep=""))
      comb.r <- c(comb.r, as.character(primer[i, "specific_reverse_primer"]))
    }
  }
  primer <- cbind(primer, f_primer = comb.f, r_primer = comb.r)
  write.csv(primer, file = snakemake@output[[1]], row.names = FALSE)
}
