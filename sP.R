####################################################################################################################################################################
############################################# ENTRY POINT - FUNCTION CALL NEW PIPELINE   ###########################################################################
####################################################################################################################################################################
command = commandArgs(trailingOnly = FALSE)
arg = commandArgs(trailingOnly = TRUE)
if(length(arg)> 0){
  conf = arg[1]
}else{
  conf = "config.conf"
}
cat(paste(Sys.time(), "| This config-file is in use: ", conf, "\n", sep =""))

source ('./bin/fun.R')

cat(paste(Sys.time(), "| Starting to fill pipeline....\n", sep =""))
nP <- newPipeline(
  command.line = command,
  ### the first entry in the first field of the config.file specifies the name of the input file
  input_folder = gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep filename ", basename(conf)), intern = TRUE)),
  type = "illumina",
  with_primertable=TRUE,
  primerTable = paste(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep filename ", basename(conf)), intern = TRUE)), ".csv", sep=""),
  ### input-file und primertable file have the same name, but primer Table is a .csv-File
  cut_off_base=as.numeric(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep basequal ", basename(conf)), intern = TRUE))),
  cut_off_seq=as.numeric(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep meanqual ", basename(conf)), intern = TRUE))),
  similarity=as.numeric(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep clustering ", basename(conf)), intern = TRUE))),
  chimera = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep chim_rm ", basename(conf)), intern = TRUE))),
  multicore = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep cores ", basename(conf)), intern = TRUE))),
  min_length_of_sequence = as.numeric(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep minlen ", basename(conf)), intern = TRUE))),
  with_tax = TRUE,
  name_of_log_file = paste("log_",format(Sys.time(), "%Y%m%d_%H%M"),".txt", sep=""), 
  name_extension = as.character(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep name_ext ", basename(conf)), intern = TRUE))),
  max.length = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep maxlen ", basename(conf)), intern = TRUE))),
  pairedEnd = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep paired_End ", basename(conf)), intern = TRUE))),
  threshold = as.double(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep threshold ", basename(conf)), intern = TRUE))),
  minoverlap = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep minoverlap ", basename(conf)), intern = TRUE))),
  minqual = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep minqual ", basename(conf)), intern = TRUE))),
  minh = as.double(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep minh ", basename(conf)), intern = TRUE))),
  mindiffs = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep mindiffs ", basename(conf)), intern = TRUE))),
  mindiv = as.double(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep mindiv ", basename(conf)), intern = TRUE))),
  beta = as.double(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep beta ", basename(conf)), intern = TRUE))),
  pseudo_count = as.double(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep pseudo_count ", basename(conf)), intern = TRUE))),
  abskew = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep abskew ", basename(conf)), intern = TRUE))),
  e_value = as.double(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep blast_Evalue ", basename(conf)), intern = TRUE))),
  filtering = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep filtering ", basename(conf)), intern = TRUE))),
  filter.method = as.character(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep filter_method ", basename(conf)), intern = TRUE))),
  filter.cutoff = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep cutoff ", basename(conf)), intern = TRUE))),
  p.correction = as.character(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep ampli_corr ", basename(conf)), intern = TRUE))),
  saving.format = as.character(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep save_format ", basename(conf)), intern = TRUE))),
  plot.ampliconduo = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep plot_AmpDuo ", basename(conf)), intern = TRUE))),
  MID.removed = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep mid_rm ", basename(conf)), intern = TRUE))),
  primer.offset = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep primer_offset ", basename(conf)), intern = TRUE))),
  ALL.removed = as.character(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep all_primer ", basename(conf)), intern = TRUE))),
  pre.filter = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep prefilter ", basename(conf)), intern = TRUE))),
  pre.filter.mq = as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep mq ", basename(conf)), intern = TRUE))),
  occurance= as.integer(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep singleton_filter ", basename(conf)), intern = TRUE))),
  cd_hit=as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep cd_hit ", basename(conf)), intern = TRUE))),
  multiqc = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep multiqc ", basename(conf)), intern = TRUE))),
  mid_check = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep mid_check ", basename(conf)), intern = TRUE))),
  skipping = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep skipping ", basename(conf)), intern = TRUE))),
  swarm = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep swarm ", basename(conf)), intern = TRUE))),
  spoint = as.logical(gsub("(.* = )(.*)( #.*$)", "\\2", system(paste("grep spoint ", basename(conf)), intern = TRUE)))
)

#### Starting point menu ####
if(nP$spoint){
  cat('Please choose the starting point of the pipeline\n1 = multiqc and primer/midcheck\n2 = assembling\n3 = dereplication\n4 = chimera_removal\n5 = filtering\n')
  interact <- readLines('stdin', n=1)
  if (interact %in% c(1:5) == F){
    stop('Only numbers from 1 to 4 are valid selections')
  }
}

#### Step 0 #####
log_folder0 <- paste0(nP$input_folder, '/logging/multiqc')

if(nP$skipping == T && file.exists(paste0(log_folder0, '/', 'multiqc_success')) == T || (nP$spoint && interact > 1)){
  cat('Status file already exists, skipping multiqc...\n')
}else{
  log_file0 <- paste0(log_folder0, '/', 'multiqc_log.txt')
  file.create(log_file0)
  
  if(nP$multiqc){
    Quality_check(nP$input_folder)
  }
  
  if(nP$mid_check){
    MID_primer_check(nP$input_folder, nP$primerTable)
  }
  file.create(paste0(log_folder0, '/', 'multiqc_success'))
  gc()
}
##############################
####   Start of Step 1    ####
##############################
## quality filtering and primer removal

##


log_folder1 <- paste0(nP$input_folder, '/logging/assembling')

if(nP$skipping == T && file.exists(paste0(log_folder1, '/', 'assembling_success')) == T || (nP$spoint && interact > 2)){
  cat('Status file already exists, skipping assembling...\n')
}else{
  log_file1 <- paste0(log_folder1, '/', 'assembling_log.txt')

  if(nP$pairedEnd){
    #### assembly of paired end reads
    #### function can be found in fun.R
    read.Assembler(nP$input_folder, nP$threshold, nP$minoverlap, nP$minqual, nP$min_length_of_sequence, nP$max.length, 
                  nP$primer.offset, nP$pre.filter, nP$pre.filter.mq)
  }else{ ## singleEnd
    source("./bin/apply.DigDeeper1.03.R")
    ## load primer table
    primer  <- read.csv(nP$primerTable, header=T, stringsAsFactors = F)
    comb = c()
    primer <- cbind(primer, combined_primer = comb)
    #folder = dirname(input_folder)
    files = list.files(input_folder, full.names=TRUE, pattern="\\.fastq$")
    
    ### apply.DigDeeper, iterate over all files
    for(i in (1:length(files))){
      file <- files[[i]]
      file.name <- unlist(strsplit(file, split = "/"))
      file.name <- file.name[length(file.name)]
      file.name <- unlist(strsplit(file.name, split = ".fastq" ))
      ## get the MID-Primer sequence for the sample
      file.name.no <- unlist(strsplit(file.name, split = name_extension ))
      this.primer <- primer[primer$Probe_FWD == file.name.no,]
      this.primer <- as.character(this.primer$combined_primer)
      
      apply.DigDeeper(file, file.name, this.primer, nP.max.length )
    }
    rm(files)
  }
  gc()
}

##############################
####    End of Step 1     ####
##############################


##############################
####   Start of Step 2    ####
##############################
## DEREPLICATION 

log_folder2 <- paste0(nP$input_folder, '/logging/dereplication')

if(nP$skipping == T && file.exists(paste0(log_folder2, '/', 'dereplication_success')) == T || (nP$spoint && interact > 3)){
  cat('Status file already exists, skipping dereplication...\n')
}else{
  log_file2 <- paste0(log_folder2, '/', 'dereplication_log.txt')
  
  ################ with cd-hit ###################
  if(nP$cd_hit==T){
    cat("cd_hit:", nP$cd_hit, "\n")
    cat(paste(Sys.time(), "| Starting with cd-hit...\n", sep=""))
    
    dirs.names <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = F)
    dirs.names <- dirs.names[! grepl('multiqc', dirs.names, fixed = T)]
    dirs <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = T)
    dirs <- dirs[! grepl('multiqc', dirs, fixed = T)]
    cat("Checking for dirs to use: \n", dirs, "\n\n")
    cat(paste(Sys.time(), "| Start dereplicating ... \n", sep=""))
    for(i in 1:length(dirs)) {
      
      f <- paste(dirs[i], "/", dirs.names[i], ".fasta", sep = "" )
      o <- sub(".fasta", ".clustered100.fasta", f)
      cdhit.log <- cdhit.clustering100(input=f, output=o, cores=nP$cores, th=1.0)
    }
    cat(paste(Sys.time(), "| END dereplicating ... \n\n", sep=""))
  }
  ################   end with cd-hit   ############
  ## read all .fasta files from result folder and dereplicate
  ## in case of pairedEnd processing, shorter reads are not sorted to longer ones
  
  if(nP$cd_hit==FALSE){
    cat("cd_hit:", nP$cd_hit, "\n")
    getlabels <- function(name, frequencies){
      labels <- sapply(X = 1:length(frequencies), FUN = function(x){
        label <- paste(name,"_", x,";size=", frequencies[x], ";", sep = "")
        return(label)
      })
    }
    
    dirs.names <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = F)
    dirs.names <- dirs.names[! grepl('multiqc', dirs.names, fixed = T)]
    dirs <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = T)
    dirs <- dirs[! grepl('multiqc', dirs, fixed = T)]
    
    ### dereplicates, removes duplicate sequences, counts occurences
    ### iterates over all files
    cat(paste(Sys.time(), "| Start dereplicating ... \n", sep=""))
    for(i in 1:length(dirs)) {
      cat(i, " \n")
      f <- paste(dirs[i], "/", dirs.names[i], ".fasta", sep = "" )
      if(file.exists(f)){
        unsorted.data <- read.fasta(f, as.string = T) #  list of vector of chars, SeqFastadna
        sequences = unlist(getSequence(unsorted.data, as.string=T)) ## list, each read as a character string
        rm(unsorted.data)
        # cat("typeof(f): ", typeof(f), "\n")
      }
      # extracts sample name
      # name = unlist(strsplit(f, split = "/"))
      # name <- name[length(name)]
      # name <- sub(".fasta", "", name)
      name <- sub(".fasta", "", basename(f))
      
      ## generate data.table for fast sorting via hash
      sequence.table <- data.table(sequences, count = rep(1, times = length(sequences)))
      cluster.table <- sequence.table[ , sum(count), by = sequences]
      setnames(cluster.table, old = "V1", new = "count" )
      ## convert clustertable back to data.frame
      
      cluster.table <- as.data.frame(cluster.table)
      text <- paste("\n", Sys.time(), "| ", dirs.names[i], " has ", length(cluster.table[,1]), " unique sequences. \n \n", sep = "")
      cat(text, file = log_file2, append = T)
      cat("cluster.table: ", typeof(cluster.table))
      
      #####################  ------ check if short sequences excist  ! not with pairedEnd !----  ########################################
      if(nP$pairedEnd == FALSE){
        max.length <- max(width(as.character(cluster.table[,1])))
        min.length <- min(width(as.character(cluster.table[,1])))
        
        ### sort short sequences to longer sequences with maximum count
        if(max.length > min.length){
          cluster.width <- width(as.character(cluster.table[,1]))
          # finds short sequences
          short.cluster.indices<- which(cluster.width < max.length)
          short.cluster<- cluster.table[short.cluster.indices,] ## data frame short sequence and occurence of the sequence
          
          write.csv(short.cluster, file = paste(.nP$resultsFolder, "/", name, "/", name, "_short_derep_cluster.csv", sep ="" ))
          text <- "shorter sequences \n"
          cat(text, file = log_file2 , append = T)
          text <- paste(short.cluster.indices)
          cat(text, file = log_file2 , append = T)
          
          ## returns a list with the indices of the short sequences that matched a long sequence
          ## writes into cluster.table new frequencies
          long.match = c()
          l.match <- foreach(x = 1:length(short.cluster[,1]), .combine='c') %do% {
            ## finds indices of matching sequences
            match <- grep(pattern = as.character(short.cluster[x, 1]), x = as.character(cluster.table[,1 ]))
            ## if more then 1 sequences matched, finds the most abundant sequence
            if(length(match) > 0){
              if(length(match) > 1){
                match.seq <- cluster.table[match,]# data frame sequence + abundance
                ## checks that sequence has max length and max abundance
                add.to <- which(width(as.character(match.seq[,1]))== max.length & match.seq[,2] == max(match.seq[,2]))
                add.to <- match[add.to]
              }else{ ## length(match == 1)
                add.to <- match
              }
              ## adds occurences of short sequence to matched long sequence
              index <- short.cluster.indices[x]
              cluster.table[add.to, 2] <- (cluster.table[add.to, 2] + cluster.table[index,  2])
              long.match <- c(long.match, index)
            }#else
            ## no match for the short sequence was found
          }
          
          ## removes short sequences that were matched
          cat(long.match)
          
          ###
          if(length(long.match) > 0){
            cluster.table <- cluster.table[- long.match, ]
          }
          
          text <- paste(length(long.match), " short sequences could be matched to long sequences \n",
                        length(cluster.table[,1]), " dereplicated cluster remain \n", sep = "")
          cat(text, file = log_file2 , append = T)
          
          ## not matched short sequences
          short.nomatch.ind <- setdiff(short.cluster.indices, long.match)
          short.nomatch.seq <- cluster.table[short.nomatch.ind,]
          write.csv(short.nomatch.seq, file = paste(nP$resultsFolder, "/", name, "/", name, "_short_derep_nomatch_cluster.csv", sep ="" ))
        }else{
          cat("no short sequences to match to longer sequences \n", file = log_file2 , append = T)
        }
      } ###################################################### end short seq check, paired End = FALSE ######################################
      
      ## arrange sequences in cluster.table in descending order
      desc.order <- order(cluster.table$count, decreasing = TRUE)
      cluster.table <- cluster.table[desc.order, ]
      
      ## sequences to upper case
      upper.seqs <- as.character(cluster.table[ ,1])
      upper.seqs <- sapply(upper.seqs, toupper)
      
      labels <- getlabels(name, cluster.table$count)
      seqs <- as.list(upper.seqs)
      
      ## save file in sample specific folder
      write.fasta(sequences = seqs , names = labels,  file.out = paste(nP$resultsFolder, "/", name, "/",  name, 
                                                                      ".clustered100.fasta", sep = ""))
      
      rm(cluster.table)
      rm(upper.seqs)
      rm(labels)
      rm(seqs)
    }
  }
  file.create(paste0(log_folder2, '/', 'dereplication_success'))
  gc()
  cat("\n")
  cat(paste(Sys.time(), "| all fasta files dereplicated (full length) \n\n", sep = ""), file = log_file2 , append = TRUE)
}
##############################
####    End of Step 2     ####
##############################



##############################
####   Start of Step 3    ####
##############################
## CHIMERA REMOVAL

log_folder3 <- paste0(nP$input_folder, '/logging/chimera_removal')

if(nP$skipping == T && file.exists(paste0(log_folder3, '/', 'chimera_removal_success')) == T || (nP$spoint && interact > 4)){
  cat('Status file already exists, skipping chimera_removal...\n')
}else{
  log_file3 <- paste0(log_folder3, '/', 'chimera_removal_log.txt')
  
  ## folder for the non chimeric data
  if(nP$similarity == 100){
    final.data.folder <- paste(nP$resultsFolder, "/finalData", sep="")
    folder.sys.call = system(paste("mkdir -p ", final.data.folder , sep=""), intern = TRUE)
    cat("CHIMERA REMOVAL \n")
    cat("final.data.folder: ", final.data.folder, "\n\n")
  }
  
  ##########
  #  reads Uchime output table and extracts positive (Y) chimeras
  #  removes Chimeras from input fasta file

  dirs.names <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = F)
  dirs.names <- dirs.names[! grepl('multiqc', dirs.names, fixed = T)]
  dirs <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = T)
  dirs <- dirs[! grepl('multiqc', dirs, fixed = T)]
  
  if(nP$chimera){
    cat(paste(Sys.time(), "| Start removing chimera: \n ", sep=""), file = log_file3 , append = TRUE )
    cmd = "./bin/usearch7 -help 2>&1"
    sys = system(cmd, intern = T)
    cat(paste0(sys, sep =" \n"), file = log_file3 , append = TRUE )
    for (n in 1:length(dirs)){
      dir <- dirs[n]
      dir.name <- dirs.names[n]
      input.file <-  paste(nP$resultsFolder, "/", dir.name, "/",  dir.name, ".clustered100.fasta", sep = "")
      uchimeout.file <- paste(nP$resultsFolder, "/", dir.name, "/",  dir.name, ".clustered100.uchime.txt", sep = "")
      if(file.exists(input.file)){
        if(nP$similarity == 100){
          nonchimera.file <- paste(final.data.folder , "/",  dir.name, ".clustered100.nonchimera.fasta", sep = "")
        }else{
          nonchimera.file <- paste(nP$resultsFolder, "/",  dir.name, "/",  dir.name, ".clustered100.nonchimera.fasta", sep = "")
        }
        
        chimera.file <- paste(nP$resultsFolder, "/", dir.name, "/",  dir.name, ".clustered100.chimera.fasta", sep = "")
        cat("Checking the variables: \n input:", input.file,"\n uchimeout:", uchimeout.file,"\n nonchimera:", nonchimera.file, 
            "\n chimera:", chimera.file, "\n\n")
        cmd <- paste("./bin/usearch7 -uchime_denovo ", input.file, " -minh ", nP$minh, " -mindiffs ", nP$mindiffs, " -mindiv ", 
                    nP$mindiv, " -xn ", nP$beta, " -dn " , nP$pseudo_count ," -abskew ", nP$abskew , " -uchimeout ", 
                    uchimeout.file, " -chimeras ", chimera.file, " 2>&1 ", sep = "")
        sys = system(cmd, intern = T, ignore.stderr = F)
        sys.length <- length(sys)
        sys <- sys[(sys.length - 1):sys.length]
        sys <-paste0(sys, collapse=" \n ")
        cat( paste("\n", Sys.time(), "| ", dir.name, ": \n", sys, sep=""), file = log_file3 , append = TRUE )
        
        tmp <- getPositiveChimeras(uchimeout.file, input.file, nonchimera.file)
        cat("\n\n")
      }
    }
  }else{
    for (n in 1:length(dirs)){
      dir <- dirs[n]
      dir.name <- dirs.names[n]
      input.file <-  paste(nP$resultsFolder, "/", dir.name, "/",  dir.name, ".clustered", similarity, ".fasta", sep = "")
      nonchimera.file <- paste(final.data.folder , "/",  dir.name, ".clustered", similarity, ".nonchimera.fasta", sep = "")
      cmd = paste("cp ", input.file, " ", nonchimera.file, sep = "")
      sys = system(cmd, intern = T, ignore.stderr = F)
    }
    cat( "chimera removal is off \n", file = log_file3 , append = TRUE )
  } 
  
  ########################################  CLUSTERING < 100 % ################################################################################
  #############################################################################################################################################
  ## if similarity < 100 -> clustering usearch deterministisch  !! don't use with  pairedEnd
  
  if(nP$similarity < 100){
    
    text <- paste(Sys.time(), " Start clustering to ", similarity, "% identity \n", sep = "")
    cat(text , file = log_file3 , append = TRUE)
    cat("start clustering: ", similarity, "% indentity \n" )
    
    dirs.names <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = F)
    dirs <- list.dirs(path = nP$resultsFolder, recursive = F, full.names = T)
    identity = as.numeric(similarity)/100
    
    final.data.folder <- paste(nP$resultsFolder, "/finalData", sep="")
    folder.sys.call = system(paste("mkdir -p ", final.data.folder , sep=""), intern = TRUE)
    
  
    ##iterates over samples
    for (k in 1:length(dirs)){
      dir <- dirs[k]
      dir.name <- dirs.names[k]
      
      text <- paste("\n", Sys.time(), " Start clustering ", dir.name , " \n", sep = "")
      
      input.file <- paste( dir, "/", dir.name , ".clustered100.nonchimera.fasta", sep = "")
      saving.file <- paste(dir, "/", dir.name, "_clustered_", identity, sep = "")
      
      ##make a new cluster folder
      cluster.folder <- paste(dir, "/cluster", sep = "")
      system(paste("mkdir -p ", cluster.folder, sep=""), intern = TRUE)
      
      ## cluster with usearch 7
      cluster.cmd = paste("./bin/usearch7 -cluster_smallmem ", input.file, " -id ", identity,
                          " -centroids ", saving.file, ".centroids",
                          " -clusters ", dir, "/cluster/", dir.name, "_clust",
                          " -uc ", saving.file, ".txt -fulldp 2>&1",  sep = "")
      
      sys = system(cluster.cmd, intern=TRUE, ignore.stderr = F)
      
      if(k == 1){
        sys1 <-paste0(sys[1:5], collapse=" \n ")
        cat(paste("usearch version: \n", sys1 , "\n", sep=""), file =  log_file3 , append = TRUE)
      }
      
      cat(text , file = log_file3 , append = TRUE)
      sys.length <- length(sys)
      sys <- sys[(sys.length - 10):sys.length]
      sys <-paste0(sys, collapse=" \n ")
      cat( sys, file = log_file3 , append = TRUE )
      
      ## iterates over all clusters, calculates size of clusters
      clusters <- list.files(cluster.folder, full.name = TRUE)
      
      ## calculates cluster size
      for( l in 1:length(clusters)){
        tmp.centroid <- count.cluster.size(clusters[l], dir, dir.name)
      }
      
      ## sort by size
      cmd <- paste("./bin/usearch7 -sortbysize ", paste(dir, "/", dir.name, "_","centroids.abundance.fasta ", sep = ""), "-output ",
                  paste(final.data.folder, "/", dir.name, ".clustered", similarity, ".nonchimera.fasta", sep = ""), " -minsize 1", sep = "")
      sys = system(cmd, intern = T)
      
    }
    
    cat(paste( "\n \n", Sys.time(),"  clustering done. \n \n", sep = ""), file = log_file3 , append = TRUE )
  }
  file.create(paste0(log_folder3, '/', 'chimera_removal_success'))
  gc()
}
##############################
####    End of Step 3     ####
##############################


##############################
####   Start of Step 4    ####
##############################
## CREATE TABLE WITH ALL SAMPLES

log_folder4 <- paste0(nP$input_folder, '/logging/filtering')

if(nP$skipping == T && file.exists(paste0(log_folder4, '/', 'filtering_success')) == T){
  cat('Status file already exists, skipping filtering...\n')
}else{
  log_file4 <- paste0(log_folder4, '/', 'filtering_log.txt')
  sys = system(paste("mkdir -p ", nP$resultsFolder, "/tables", sep=""), intern = TRUE)
  final.data.folder <- paste(nP$resultsFolder, "/finalData", sep="")

  table_1 <- make.table(final.data.folder,  nP$input_folder)
  table_1
  
  ####################### FILTERING ################################################################################################################
  ##################################################################################################################################################
  
  cat("Start with filtering: ",nP$filtering, "\n")
  if(nP$filtering){
    cat("filter-method: ", nP$filter.method, "\n")
    ## AmpliconDuo based filtering
    ## requires always A-B, A-B, A-B
    if(nP$filter.method == "aduo"){
      #ampliconDuo.filter(table, final.data.folder, p.correction, saving.format, plot.ampliconduo )
      mixed.filtering(table_1, final.data.folder, nP$filter.cutoff, nP$p.correction, nP$saving.format, nP)
    }
    
    if(nP$filter.method =="singelton"){
      # cat("table: ", exists('table'), "\n")
      # cat("final.data.folder: ", exists('final.data.folder'), "\n")
      # cat("nP$filter.cutoff: ", exists('nP$filter.cutoff'), "\n")
      # cat("nP$occurance: ", exists('nP$occurance'), "\n")
      singelton.filter(table_1, final.data.folder, nP$filter.cutoff, nP$occurance)
    }
    
    ## mixed filtering, if not all samples have A and B
    if(nP$filter.method == "mixed"){
      mixed.filtering(table_1, final.data.folder, nP$filter.cutoff, nP$p.correction, nP$saving.format, nP)
    }
  }

  ################################################ SWARM #############################################################################
  #############################################################################################################################################
  
  if(nP$swarm){
    table.name <- paste0(final.data.folder,"/", nP$input_folder, "_Table.csv")
    fseq <- read.delim(table.name, stringsAsFactors = FALSE, sep= " ")
    merged_data <- paste0(final.data.folder, '/merged.fasta')
    # Create a fasta file from the table of the merged probes
    for(i in 1:nrow(fseq)){
      cat(sprintf('>%i;size=%i;\n%s\n', i, fseq[i,2]+fseq[i,3], fseq[i,1]), file = merged_data, append = TRUE)
    }
    swarm(final.data.folder, nP$cores)
  }

  ##############################
  ####    End of Step 4     ####
  ##############################
  file.create(paste0(log_folder4, '/', 'filtering_success'))
  cat(paste(Sys.time(), "| Processing finished. End of pipeline....\n", sep=""), file = log_file4, append = TRUE)

  # Write down how many sequences were discarded/kept per step
  l_folder <- paste0(nP$input_folder, '/logging/')
  l_file <- paste0(l_folder, nP$name_of_log_file)
  prin_seq <- system(paste0("grep 'Bad sequences' ", l_folder, "assembling/assembling_log.txt | cut -d ' ' -f5"), intern = TRUE)
  cat(sprintf('%i sequences were removed by prinseq\n', sum(as.numeric(gsub(',', '', prin_seq)))), file = l_file, append = TRUE)

  pan_seq <- system(paste0("grep NOALGN ", l_folder, "assembling/pandaseq_log.txt | cut -f4"), intern = TRUE)
  cat(sprintf('Pandaseq could not assemble %i sequences\n', sum(as.numeric(pan_seq))), file = l_file, append = TRUE)

  uchim <- system(paste0("egrep 'Writing [[:digit:]]* chimeras' ", l_folder, "chimera_removal/chimera_removal_log.txt | cut -d ' ' -f9"), intern = TRUE)
  cat(sprintf('Uchime found %i chimeras\n', sum(as.numeric(uchim))), file = l_file, append = TRUE)

  a_duo <- system(paste0("grep 'kept' ", l_folder, "filtering/filtering_log.txt | cut -d ' ' -f2"), intern = TRUE)
  cat(sprintf('Ampliconduo kept %s sequences\n', unique(a_duo)), file = l_file, append = TRUE)
  
  # Merge logfiles
  merge_logs <- paste0('cat ', nP$input_folder, '/logging/multiqc/multiqc_log.txt ', nP$input_folder, '/logging/assembling/assembling_log.txt ',
    nP$input_folder, '/logging/dereplication/dereplication_log.txt ', nP$input_folder, '/logging/chimera_removal/chimera_removal_log.txt ', nP$input_folder, '/logging/filtering/filtering_log.txt ',
    '>> ', nP$input_folder, '/logging/', nP$name_of_log_file)
  system(merge_logs)
}