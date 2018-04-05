#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| make.Table.fromSilva.R                              |\n")
#cat("|                                                     |\n")
#cat("| (c) by Anja Lange, 2014                             |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("| anja.lange@uni-due.de                               |\n")
#cat("|                                                     |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")


### make table with data from silva db, contain taxonomic information
final.table.silva <- function(input1, input2, output, found.bug){
  
  ### taxonomy has always one more line then tab has rows, first Line contains comlumn names
  col.classes <- "character"
  taxonomy <- fread(input1, header = T, colClasses=col.classes)
  taxonomy <- as.data.frame(taxonomy)
  #taxonomy = read.table(input1, header = T, fill= T, strip.white=TRUE) ## every row of the table as a character
  tab = read.table(input2, header=T)
  taxonomy.levels <- length(names(taxonomy)) -13
  
  ###############  if blast table is different from _Table.csv #################
  unblasted.seqs = FALSE
  if(found.bug){
    ## qids correspond to row of sequence, remove sequences that gave no blast results
    ## read
    file = file(paste(unlist(strsplit(input1, split = "_Table-taxonomy.csv")), ".noBlast", sep = ""))
    not.blasted <- readLines(con = file, warn = FALSE)
    close(file)
    not.blasted <- as.integer(not.blasted)
    tab.notBlasted <- tab[not.blasted,]
    ## writes table with sequences that gave no blast result
    write.table(tab.notBlasted, row.names = F, file = paste(unlist(strsplit(input1, split = "-taxonomy.csv")), ".noBlast", sep = "") )
    if(!(length(not.blasted) == 1 && not.blasted[1] == 0)){
      ## removes sequences with no blast result from tab
      tab = tab[-not.blasted, ]
      unblasted.seqs = TRUE
    }
  }
  
  #################################
  ## extracts informations for thr final table
  tab.res <- taxonomy[, c(2, 4, 3, 11 , 13:length(taxonomy[1, ]))]
  tab.res = cbind(sum = NA, tab.res)
  
  tab.res = cbind(tab, tab.res)
  
  if(unblasted.seqs){
    ### adds not blasted sequences to the end of the table
    #no.blast.count <- length(tab.notBlasted[,1])
    tab.res.notblasted = cbind(tab.notBlasted, sum = NA)
    tab.res.notblasted = cbind(tab.res.notblasted, subject.id = NA)
    tab.res.notblasted = cbind(tab.res.notblasted, alignment.length = NA)
    tab.res.notblasted = cbind(tab.res.notblasted, identity = NA)
    tab.res.notblasted = cbind(tab.res.notblasted, evalue = NA)
    tab.res.notblasted = cbind(tab.res.notblasted, stitle = NA)
    
    tmp <- foreach( i = 1:taxonomy.levels)%do%{
      tab.res.notblasted <-  cbind(tab.res.notblasted, NA)
    }
    names(tab.res.notblasted) = names(tab.res)
    
    tab.res <- rbind(tab.res, tab.res.notblasted)
  }
  
  ### calculates sums ## to do make faster with data.table
  sums <- foreach(i = 1:length(tab.res[,1]), .combine = "c") %dopar% {
    tmp = sum(tab.res[i,2:(ncol(tab))])
    cat(paste("sum", i, "/", length(tab.res[,1]), "\n"))
    return(tmp)
  }
  tab.res[,"sum"] <- sums
  
  ## determines sequence length
  sequence.length <- sapply(as.character(tab.res[,1]), nchar )
  tab.res <- cbind(sequence.length, tab.res)
  
  write.table(tab.res, file=output, quote=F, row.names=F, sep="\t")
  cat("Done.\n")
  return(tab.res)
}

#############################################################################################################################################
#######################################################  make.Tables.fromSilva  #############################################################
#############################################################################################################################################

make.Tables.fromSilva <- function(input, output, sequence.count, silva.version){
  ## using fread helps overcoming problem if taxonomy contains escape characters or others that are not correcly interpretet with read.table
  ## Warnmeldung:
  ## In scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  :
  ##            EOF within quoted string
  
  table.res = fread(input, header=F) ## table with blast results
  table.res = as.data.frame(table.res) ## 
  colnames(table.res) = c("query.id", "subject.id", "identity", "alignment.length", "mismatches",
                          "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score", "stitle")
  
  #######################################  bug fix and no good match found #########################################################
  ### test if blast returned same number of sequences that were queried
  cat("checking if blast was affected by blastn bug (blast version >2.2.27) \n", file = .GlobalEnv[["logfile"]] , append = TRUE)
  blast.bug <- FALSE
  ref.seq.names <- seq(1:sequence.count)
  blast.qids <- table.res[,1]
  if(sequence.count != length(blast.qids)){
    cat("sequence count and blast result count are different! \n ", file = .GlobalEnv[["logfile"]] , append = TRUE)
    blast.bug = TRUE
  }else{
    ## test if qids match
    if(!identical(ref.seq.names, blast.qids)){ ##TODO maybe better check if each element is 0
      cat("order of sequences differs between sequences and blast \n", file = .GlobalEnv[["logfile"]] , append = TRUE)
      blast.bug = TRUE
    }
  }
  
  if(blast.bug == TRUE){
    ## identify missing sequences (not a bug, if e-value was too low)
    missing.qids <- setdiff(ref.seq.names, blast.qids)
    if(length(missing.qids) > 0){## some sequences gave no blast result
      ## remove later from _table.csv table and write to seperate file
      cat("no blast result for some sequences! \n", file = .GlobalEnv[["logfile"]] , append = TRUE)
      file = paste(unlist(strsplit(input, split = "_Table.blast")), ".noBlast", sep = "")
      ids.to.file <- paste0(missing.qids, collapse= "\n")
      cat(ids.to.file, file = file)
    }else{
      file = paste(unlist(strsplit(input, split = "_Table.blast")), ".noBlast", sep = "")
      cat("0", file = file)
    }
    ### identify duplicated qids
    duplicated.qids <- duplicated(blast.qids) 
    duplicated.qids <- which(duplicated.qids == TRUE)
    if(length(duplicated.qids) > 0){
      cat("some qids were duplicated in the blast result \n", file = .GlobalEnv[["logfile"]] , append = TRUE)
      file = paste(unlist(strsplit(input, split = "_Table.blast")), ".dupBlast", sep = "")
      ids.to.file <- paste0(duplicated.qids, collapse= "\n")
      cat(ids.to.file, file = file)
      
      ## removes duplicated blast results from table
      table.res <- table.res[-duplicated.qids, ]
    }    
  }else{
    cat("not affected by blast bug! or missing blast results \n", file = .GlobalEnv[["logfile"]] , append = TRUE)
  }

  ################################### end bug fix ##########################################################
  
  ## extract taxonomy
  if(silva.version== "115"){
    taxonomy.table <- resolve.silva.taxonomy.mc(as.character(table.res$stitle))
  }else{
    taxonomy.table <- resolve.silva.taxonomy.mc.v119(as.character(table.res$stitle))
  }
  
  
  ## removes empty taxonomic ranks
  Na.vector <- as.character(rep(NA, times = length(taxonomy.table[,1])))
  empty = c()
  em <- foreach(x = 1:(length(taxonomy.table[1,]))) %do% {
    tmp <- as.character(taxonomy.table[, x])
    tmp <- identical(tmp, Na.vector)
    if(tmp == TRUE){
      empty= c(empty, x)
    }
  }
  
  if(length(empty) > 0){
    taxonomy.table <- taxonomy.table[, - empty]
  }
 
  ## removes column with colon seperated taxonomy
  #table.res <- table.res[, - 13]
  
  ### adds columns with seperated taxonomic rankes
  table.res = cbind(table.res, taxonomy.table)
  
  write.table(table.res, file=output, quote=F, row.names=F, sep="\t")
  
  ### information is required for the final.table function call
  return(blast.bug)  
}

###################################################################################################################################
########################################### resolve.silva.taxonomy - sequential version ############################################
#tax_ranks_ssu_version.csv
#-------------------------
# This file contains taxonomic rank designations for all taxonomic paths
#used in the silva ssu taxonomy.

#The file is tab delimited and uses unix line ends. There are four
#fields:

#  path:
#  The full taxonomic path including the name of the group itself. i
#Segments are separated with ";"
#node:
#  The name of the group (repeated from path)
#rank:
#  The rank designation. 
#Can be: domain, phylum, class, order, family, genus
#remark:
#  Can be empty ('') or a or w.
#a: Marks taxa of environmental origin. That is, taxa containing no 
#sequence coming from a cultivated organism.
#w: Marks taxa scheduled for revision in the next release.

resolve.silva.taxonomy <- function(taxa){ 
  require(data.table)
  tax.ranks <- fread(.GlobalEnv[["silva.tax.ranks"]], header = T)
  setkey(tax.ranks , key= "node")
  levels <- tax.ranks[, sum(), by = rank]
  levels <- levels[,rank] 
  ## all so far possible taxonomix levels
  all.ranks <- c( "domain", "superkingdom", "major_clade", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass",
                  "superorder", "order","suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tibe", "subtribe", "genus", 
                  "subgenus", "species group" , "species subgroup", "species", "subspecies", "varietas", "forma" )
  
  ordered.ranks <- intersect(all.ranks, levels)
  ordered.ranks <- c(ordered.ranks, "species")
  
  taxonomy <- as.data.frame(matrix(NA, nrow = length(taxa), ncol =length(ordered.ranks), dimnames = list(c(1:length(taxa)), ordered.ranks)))
  
  taxo <- foreach(i = 1:length(taxa)) %do% {
    tmp <- unlist(strsplit(taxa[i], split = ";"))
    foreach(j = 1:length(tmp)) %do% {
      rank <- tax.ranks[tmp[j]]$rank
      if(tmp[j] != "uncultured" && length(rank) == 1 && !is.na(rank[1])){
        taxonomy[i,rank] <- tmp[j]        
      }else{
        ## rank is not distinct, needs to consider whole path
        if(!is.na(rank[1])){
          p = paste("^", paste0(tmp[1:j], collapse=";"), "$", sep ="")
          index = grep(p, tax.ranks$path)
          rank <- tax.ranks[index, ]$rank
          taxonomy[i,rank] <- tmp[j]  
        }else{
          ## species are not in tax_ranks_ssu_version.csv
          ## if last rank is not found in table, it is assumed to be species
          if(j == length(tmp)){
            rank = "species"
            taxonomy[i,rank] <- tmp[j]
          } 
        }
      }
    }
  } 
  return(taxonomy)
}

######################################################################################################################################
########################################### resolve.silva.taxonomy.mc - parallel version ############################################
## compared parallel and sequential methods with 4 cores , taxa with 575 entries
## sequential
## User      System verstrichen 
## 17.757       0.000      17.790 
## 4 cores
## User      System verstrichen 
## 18.485       0.240       5.254 
## 8 cores
## User      System verstrichen 
##  45.451       0.476       5.096  

resolve.silva.taxonomy.mc <- function(taxa){
  require(data.table)
  tax.ranks <- fread(.GlobalEnv[["silva.tax.ranks"]], header = T)
  setkey(tax.ranks , key= "node")
  
  levels <- tax.ranks[, sum(), by = rank]
  levels <- levels[,rank]
  
  ## all so far possible taxonomix levels
  all.ranks <- c( "domain", "superkingdom", "major_clade", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass",
                  "superorder", "order","suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tibe", "subtribe", "genus", 
                  "subgenus", "species group" , "species subgroup", "species", "subspecies", "varietas", "forma" )
  
  ordered.ranks <- intersect(all.ranks, levels)
  ordered.ranks <- c(ordered.ranks, "species")
  
  taxo <- foreach(i = 1:length(taxa), .combine = 'rbind') %dopar% {
    tmp <- unlist(strsplit(taxa[i], split = ";"))
    tmp.tax <- as.data.frame(matrix(NA, nrow = 1, ncol =length(ordered.ranks), dimnames = list(i, ordered.ranks)))
    foreach(j = 1:length(tmp)) %do% {
      rank <- tax.ranks[tmp[j]]$rank
      if(tmp[j] != "uncultured" && length(rank) == 1 && !is.na(rank[1])){
        tmp.tax[1,rank] <- tmp[j]     
      }else{        
        if(!is.na(rank[1])){
          p = paste("^", paste0(tmp[1:j], collapse=";"), "$", sep ="")
          index = grep(p, tax.ranks$path)
          rank <- tax.ranks[index, ]$rank
          tmp.tax[1,rank] <- tmp[j]  
        } else{
          ## species are not in tax_ranks_ssu_version.csv
          ## if last rank is not found in table, it is assumed to be species
          if(j == length(tmp)){
            rank = "species"
            tmp.tax[1,rank] <- tmp[j]
          }            
        }
      }   
    }
    value <- tmp.tax 
  }
  return(taxo)
}
#######################################################################################################################################
##################################################### resolve.silva.taxonomy.mc.v119  #################################################
#######################################################################################################################################

resolve.silva.taxonomy.mc.v119 <- function(taxa){
  require(data.table)
  tax.ranks <- fread(.GlobalEnv[["silva.tax.ranks"]], header = F)
  setnames(tax.ranks, old = c("V1", "V2", "V3", "V4", "V5"), new = c("path", "id", "rank", "remark", "release"))

  norank <- which(tax.ranks$rank == "")  
  tax.ranks <- tax.ranks[norank, rank:= "not_specified"]
  
  ## add node to table
  nodes <- foreach(k = 1: length(tax.ranks$path), .combine = 'c') %dopar% {
    tmp = unlist(strsplit(tax.ranks$path[k], split = ";"))
    tmp = tmp[length(tmp)]
  }
  tax.ranks <- cbind(tax.ranks, node = nodes)
  
  setkey(tax.ranks , key= "node")
  
  levels <- tax.ranks[, sum(), by = rank]
  levels <- levels[,rank]
  
  ## all so far possible taxonomix levels
  all.ranks <- c( "domain", "superkingdom", "major_clade", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass",
                  "superorder", "order","suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tibe", "subtribe", "genus", 
                  "subgenus", "species group" , "species subgroup", "species", "subspecies", "varietas", "forma" , "not_specified")
  
  ordered.ranks <- intersect(all.ranks, levels)
  ordered.ranks <- c(ordered.ranks, "species")
  
 
  taxo <- foreach(i = 1: length(taxa), .combine = 'rbind') %dopar% {
    tmp <- unlist(strsplit(taxa[i], split = ";"))
    tmp.tax <- as.data.frame(matrix(NA, nrow = 1, ncol =length(ordered.ranks), dimnames = list(i, ordered.ranks)))
    x <- foreach(j = 1:length(tmp)) %do% {
      rank <- tax.ranks[tmp[j]]$rank
      if(tmp[j] != "uncultured" && length(rank) == 1 && !is.na(rank[1])){
        tmp.tax[1,rank] <- tmp[j]     
      }else{        
        if(!is.na(rank[1])){
          p = paste("^", paste0(tmp[1:j], collapse=";"), "$", sep ="")
          index = grep(p, tax.ranks$path)
          rank <- tax.ranks[index, ]$rank
          tmp.tax[1,rank] <- tmp[j]  
        } else{
          ## species are not in tax_ranks_ssu_version.csv
          ## if last rank is not found in table, it is assumed to be species
          if(j == length(tmp)){
            rank = "species"
            tmp.tax[1,rank] <- tmp[j]
          }            
        }
      }   
    }
    value <- tmp.tax 
  } 
  return(taxo)
}
