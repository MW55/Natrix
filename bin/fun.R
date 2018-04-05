#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| newPipeline_main_v0.5.R                             |\n")
#cat("|                                                     |\n")
#cat("| (c) by Anja Lange and Dominik Heider, 2012-2014     |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("| anja.lange@uni-due.de                               |\n")
#cat("|                                                     |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")

### New pipeline, for sorted files
#if(require("") == F) {
#  install.packages("") 
#  require("")
#}
## check wich packages are from bioconductor
require(Rsamtools)
require(GenomicRanges)
require(ShortRead)
require(Biostrings)
require(IRanges)
require(BiocGenerics)
require(parallel)
require(doMC)
require(caTools)
require(seqinr)
require(data.table)
require(utils)

newPipeline <- function(
  command.line,
  input_folder,
  ### the path to the input folder, no slash at the end
  type="illumina",
  ### the type of input data. Currently illumina and 454
  with_primertable=FALSE,
  ### if there is a primertable for sorting. If the samples are already sorted, then there is probably none
  primerTable = "",
  ### the table with the primers, including path if not in wd
  cut_off_base=15,
  ### the cut off for a single base to kill the sequence
  cut_off_seq=25,
  ### the cut off for the whole sequence to kill the sequence (mean)
  similarity,
  ### the similarity to put sequences together inside a cluster
  chimera,
  ### if chimeras should be removed
  multicores,
  ### number of cores to use, if not specified will take all cores it detects
  min_length_of_sequence = 0,
  ### the minimum length a sequence should have top be processed further
  with_tax = TRUE,
  ### if taxonomy search should be used. Only in combination with blast!
  name_of_log_file = "logging.txt",
  ###
  name_extension = "",
  ## difference between Probe_FWD and file name
  max.length,
  ##
  pairedEnd,
  ## specifies if forward and reverse reads are supplied
  threshold,
  ## only used if pairedEnd = TRUE, score treshold used for pandaseq
  minoverlap,
  ## only used if pairedEnd = TRUE, minimal overlap between forward and reverse reads
  minqual,
  ## only used if pairedEnd = TRUE, minimal quality value for bases in a assembled read to be accepted
  minh,
  ## minimum score for uchime_denovo
  mindiffs,
  ## minimum number of diffs in a segment
  mindiv,
  ## minimum divergence
  beta,
  ##weight of no vote
  pseudo_count,
  ##Pseudocount prior for no votes
  abskew,
  ## multiplicity of parent over chimera
  filtering,
  ## decides wether filtering should be applied
  filter.method,
  ## specifies the filtering method, ampduo or cutoff
  filter.cutoff,
  ## cutoff for cutoff filtering
  p.correction,
  ## for ampliconduo
  saving.format,
  ## for ampliconduo.plot
  plot.ampliconduo,
  ##
  MID.removed,
  ## ncbi-blast-2.2.29+
  primer.offset,
  ##
  ALL.removed,
  ## only with pandaseq
  pre.filter, 
  ## mean quality treshold
  pre.filter.mq,
  ## continues pipeline with the blast
  occurance,
  # logical, indicates if cd-hit should be used for dereplication
  cd_hit, 
  ##
  multiqc,
  ##
  mid_check,
  ##
  skipping,
  ##
  swarm,
  ##
  spoint,
  ...
){

  ##### function returns a character vector with all directories
  cmd = paste("find " ,input_folder, "/ -name *.gz", sep = "")
  if(system(cmd)){
    cat(paste(Sys.time(), "| no *.gz file found => nothing to unzip \n", sep=""))
  }else{
    cat(paste(Sys.time(), "| found a gz-file... \n", sep=""))
    unzip.cmd <- paste( "for i in ", input_folder,"/*.gz; do (gzip -kd $i) done", sep = "")
    unzip.sys.call = system(unzip.cmd, intern=TRUE)
  }
   

  #### preparations for parallel computing
    if(is.na(multicores)){
    cores <- detectCores()
  }else{
    cores <- multicores
  }
  cat(paste(Sys.time(), "| check for the use of multiple cores => ", cores, "\n", sep=""))

  os <- .Platform$OS.type
  if (os == "unix"){
    require(parallel)
    require(doMC)
    registerDoMC(cores)
    ### registerDoMC is required to execute foreach in parallel, if only the doMC package is loaded, but the
    ### the registerDoMC method is not called, the foreach loop will run sequential
  }else{
    require (snow)
    require (doSNOW)
    cl <- makeCluster(cores, type = "SOCK", outfile = "")
    registerDoSNOW(cl)
    clusterEvalQ(cl, c(require(ShortRead), require(seqinr), require(parallel), require(Biostrings), require(foreach)))
    clusterExport(cl, c("qualy_apply", "cluster_apply", "cluster_apply_inner"))
  }
  cat(paste(Sys.time(), "| check of Platform type => ", os, "\n", sep=""))
  
  resultsFolder <- paste(input_folder, "/results", sep = "" )
  
  ### logging und result-directories are created
  ### system(command, ..) invokes the OS command specified by command
  ### paste Concatenates vectors after converting to character, sep defines the seperation of the terms
  sys = system(paste("mkdir -p ",input_folder, "/logging", sep=""), intern = TRUE)
  sys = system(paste("mkdir -p ",input_folder, "/logging/multiqc", sep=""), intern = TRUE) 
  sys = system(paste("mkdir -p ",input_folder, "/logging/assembling", sep=""), intern = TRUE)
  sys = system(paste("mkdir -p ",input_folder, "/logging/dereplication", sep=""), intern = TRUE)
  sys = system(paste("mkdir -p ",input_folder, "/logging/chimera_removal", sep=""), intern = TRUE)
  sys = system(paste("mkdir -p ",input_folder, "/logging/filtering", sep=""), intern = TRUE)
  sys = system(paste("mkdir -p ",resultsFolder, sep=""), intern = TRUE)
  sys = system(paste("mkdir -p ",resultsFolder, "/multiqc", sep=""), intern = TRUE)

  if(is.null(input_folder)==FALSE){ ### tests for the existence of an input-file
    cat(paste("\n", Sys.time(), "| Input folder '", input_folder, "' is existing \n", sep=""))
    log_file <<- paste(input_folder, "/logging/", name_of_log_file, sep="")
    if(file.exists(log_file)==FALSE){
      cat(paste(Sys.time(), "| There is no logfile with name '", log_file, "'; creating it ...\n\n", sep=""))
      file.create(log_file)
    }else{
      cat(paste(Sys.time(), "| Logfile '", log_file, "' is already existing. \n Saving that logfile to new file with ending '_old'. \n", sep=""))
      file.copy(log_file, paste(log_file, "_old", sep=""), overwrite = TRUE, copy.mode = TRUE)
      cat(paste(Sys.time(), "| New logfile '", log_file, "' will be created. \n\n", sep=""))
      file.create(log_file)
    }
    ##### START: to check for a specific logfile #####
    # log_file_test <<- "./test/logging/log_20161220.txt"  
    # if(file.exists(log_file_test)==FALSE){
    #   cat(paste(Sys.time(), "| There is no log file with name", log_file_test, "...\n", sep=""))
    # }else{
    #   cat(paste(Sys.time(), "| logfile", log_file_test, "is already existing. \n", sep=""))
    # }
    ##### END: to check for a specific logfile #####
 }else{
    cat(paste(Sys.time(), "| Input_Folder ", input_folder, " is not existing \n\n", sep=""))
 }

  cat(paste("\n", Sys.time(), "| These parameters were provided to pipeline: \n", sep=""), file = log_file, append = T)
  cat(command.line, "\n", file = log_file , append = T)
  cat("Input_folder: ", input_folder, "\n", file = log_file , append = T)
  ## writes configuration file parameter into logfile
  cat(paste(" \n number of cores: ", cores), file = log_file , append = T)
  cat(paste(" \n paired end assembly:", pairedEnd), file = log_file , append = T)
  cat(paste(" \n minlength: ", min_length_of_sequence), file = log_file , append = T)
  cat(paste(" \n max_length: ", max.length), file = log_file , append = T)

  if(pairedEnd){
    cat(paste(" \n\n prefiltering:", pre.filter ), file = log_file, append = T)
    if(pre.filter){
      cat(paste(" \n minimum mean quality: ", pre.filter.mq), file = log_file, append = T)
    }
    cat(paste(" \n\n parameter for PandaSeq: "), file = log_file, append = T)
    cat(paste(" \n threshold: ", threshold), file = log_file, append = T)
    cat(paste(" \n minoverlap: ", minoverlap), file = log_file, append = T)
    cat(paste(" \n minqual: ", minqual), file = log_file, append = T)
    cat(paste(" \n primer offset: ", primer.offset), file = log_file, append = T)
    cat(paste(" \n MID removed: ", MID.removed), file = log_file, append = T)
    cat(paste(" \n ALL removed: ", ALL.removed), file = log_file, append = T)
   
    
  }else{
    cat(paste(" \n basequality: ", cut_off_base), file = log_file, append = T)
    cat(paste(" \n meanquality: ", cut_off_seq, "\n \n"), file = log_file, append = T)
  }

  cat(paste("\n\n cd_Hit clustering = ", cd_hit), file = log_file, append = T)
  
  if(chimera){
    cat(paste("\n\n parameter for chimera removal: "), file = log_file , append = T)
    cat(paste(" \n minh: ", minh), file = log_file , append = T)
    cat(paste(" \n mindiffs: ", mindiffs), file = log_file , append = T)
    cat(paste(" \n mindiv: ", mindiv), file = log_file , append = T)
    cat(paste(" \n abskew: ", abskew), file = log_file , append = T)
    cat(paste(" \n pseudo_count: ", pseudo_count), file = log_file , append = T)
    cat(paste(" \n beta: ", beta), file = log_file , append = T)
    cat("\n", file = log_file , append = T)
  }else{
    cat(paste(" \n Uchime: OFF "), file = log_file , append = T)
  }

  if(filtering){
    cat(paste(" \n filtering : ", filter.method), file = log_file , append = T)
    if(filter.method != "aduo"){
      cat(paste(" \n cutoff: ", filter.cutoff), file = log_file , append = T)
    }
  }else{
    cat(paste(" \n Filter: OFF "), file = log_file , append = T)
  }

  cat(paste("\n \n clustering: ", similarity, "\n \n"), file = log_file , append = T)

  
  ##### to make all variables available in sP.R #####
  nP <- list(command.line = command.line,
             input_folder = input_folder,
             type = type,
             with_primertable=with_primertable,
             primerTable=primerTable,
             cut_off_base=cut_off_base,
             cut_off_seq=cut_off_seq,
             similarity=similarity,
             chimera=chimera,
             multicores=multicores,
             min_length_of_sequence=min_length_of_sequence,
             with_tax=with_tax,
             name_of_log_file=name_of_log_file,
             name_extension=name_extension,
             max.length=max.length,
             pairedEnd=pairedEnd,
             threshold=threshold,
             minoverlap=minoverlap,
             minqual=minqual,
             minh=minh,
             mindiffs=mindiffs,
             mindiv=mindiv,
             beta=beta,
             pseudo_count=pseudo_count,
             abskew=abskew,
             filtering=filtering,
             filter.method=filter.method,
             filter.cutoff=filter.cutoff,
             p.correction=p.correction,
             saving.format=saving.format,
             plot.ampliconduo=plot.ampliconduo,
             MID.removed=MID.removed,
             primer.offset=primer.offset,
             ALL.removed=ALL.removed,
             pre.filter=pre.filter,
             pre.filter.mq=pre.filter.mq,
             occurance=occurance,
             cd_hit=cd_hit,
             log_file=log_file,
             cores=cores,
             os=os,
             resultsFolder = resultsFolder,
             multiqc = multiqc,
             mid_check = mid_check,
             skipping = skipping,
             swarm=swarm,
             spoint=spoint)
  return(nP)
}  

####### Step 0 (optional) ########
Quality_check <- function(input_folder){
  # Initial quality check with fastqc & multiqc
  multiqc <- "./bin/MultiQC"
  fastqc.cmd <- paste0('find ', input_folder, '/ -name *.gz', ' | ', 'xargs ', 'fastqc -o ', 
    nP$resultsFolder, '/multiqc')
  system(fastqc.cmd)
  multiqc.cmd <- paste0('multiqc . -o ', nP$resultsFolder, '/multiqc')
  system(multiqc.cmd)
}

MID_primer_check <- function(input_folder, primerTable){
  # Check for primers & MIDs in the sequences
  primer <- read.csv(nP$primerTable, header=T, stringsAsFactors = F)
  files_fwd <- list.files(input_folder, full.names=TRUE, pattern="\\R1.fastq.gz$")
  files_rev <- list.files(input_folder, full.names=TRUE, pattern="\\R2.fastq.gz$")

  # Change the IUPAC nucleotide code to corresponding regular expressions
  regex_sub <- function(input_sequence){
    input_sequence <- gsub('Y', '[CT]', input_sequence)
    input_sequence <- gsub('R', '[AG]', input_sequence)
    input_sequence <- gsub('S', '[GC]', input_sequence)
    input_sequence <- gsub('W', '[AT]', input_sequence)
    input_sequence <- gsub('K', '[GT]', input_sequence)
    input_sequence <- gsub('M', '[AC]', input_sequence)
    input_sequence <- gsub('B', '[CGT]', input_sequence)
    input_sequence <- gsub('D', '[AGT]', input_sequence)
    input_sequence <- gsub('H', '[ACT]', input_sequence)
    input_sequence <- gsub('V', '[ACG]', input_sequence)
    input_sequence <- gsub('N', '[ATCG]', input_sequence)
    return(input_sequence)
  }
  
  # Change the sequence to its reverse complementary sequence and change the IUPAC nucleotide code
  # to the corresponding COMPLEMENTARY code
  complement_seq <- function(input_sequence){
    input_sequence <- reverse(input_sequence)
    input_sequence <- chartr("ATGC", "TACG", input_sequence)
    input_sequence <- gsub('Y', '[GA]', input_sequence)
    input_sequence <- gsub('R', '[TC]', input_sequence)
    input_sequence <- gsub('S', '[CG]', input_sequence)
    input_sequence <- gsub('W', '[TA]', input_sequence)
    input_sequence <- gsub('K', '[CA]', input_sequence)
    input_sequence <- gsub('M', '[TG]', input_sequence)
    input_sequence <- gsub('B', '[GCA]', input_sequence)
    input_sequence <- gsub('D', '[TCA]', input_sequence)
    input_sequence <- gsub('H', '[TGA]', input_sequence)
    input_sequence <- gsub('V', '[TGC]', input_sequence)
    input_sequence <- gsub('N', '[ATCG]', input_sequence)
    chartr("ATGC", "TACG", input_sequence)
    return(input_sequence)
  }
  
  # Check if the forward primer and/or the MID are still in the fastq sequences
  ignore = FALSE
  for (i in 1:nrow(primer)){
    seq_sa <- readFastq(files_fwd[grep(primer[i,1], files_fwd)])
    mid_fwd <- paste0(regex_sub(primer[i,7]), regex_sub(primer[i,8]))
    len_mid_primer <- nchar(primer[i,8]) + nchar(primer[i,7])
    if(any(lapply(mid_fwd, grepl, sread(seq_sa))[[1]]) || any(lapply(regex_sub(primer[i,7]), 
      grepl, sread(seq_sa))[[1]]) || any(lapply(regex_sub(primer[i,8]), 
      grepl, sread(seq_sa))[[1]])){
        for (j in 1:100){
          seq_start <- sread(seq_sa)[[j]][1:len_mid_primer]
          if(grepl(mid_fwd, seq_start) || grepl(regex_sub(primer[i,7]), seq_start) || grepl(regex_sub(primer[i,7]), seq_start)){
            text <- sprintf('Found MID/Primer of probe %s in read %i \n\n', primer[i,1], j)
            cat(text, file=log_file0, append = TRUE)
            if (ignore == FALSE){
              cat(text)
              if(nP$ALL.removed){
                cat('Do you want to abort the process?\nPlease enter a number from the selection\n1 = Yes\n2 = No\n3 = Ignore further warnings\n')
                interact <- readLines('stdin', n=1);
                if(interact == 1){
                  stop('Pipeline stopped')
                } else if(interact == 2){
                  cat('Found primers/MIDs were written to config file')
                } else {
                ignore <- TRUE
                }
              }
            }
          }        
        } 
    }    
  }

  # Check if the reverse primer and/or the MID is still in the fastq sequences
  ignore = FALSE
  for (i in 1:nrow(primer)){
    seq_sa <- readFastq(files_rev[grep(primer[i,1], files_rev)])
    rev_mid <- paste0(complement_seq(primer[i,10]), complement_seq(primer[i,7]))
    len_mid_primer <- nchar(primer[i,7]) + nchar(primer[i,10])
    if(any(lapply(rev_mid, grepl, sread(seq_sa))[[1]]) || any(lapply(primer[i,10], 
      grepl, sread(seq_sa))[[1]]) || any(lapply(complement_seq(primer[i,10]), grepl, 
      sread(seq_sa))[[1]]) || any(lapply(complement_seq(primer[i,7]), grepl, 
      sread(seq_sa))[[1]])){
        for (j in 1:100){
          seq_end <- tail(sread(seq_sa)[[j]], n=len_mid_primer)
          if(grepl(rev_mid, seq_end) || grepl(primer[i,7], seq_end) || grepl(primer[i,10], seq_end)){
            text <- sprintf('Found MID/Primer of probe %s in read %i \n\n', primer[i,1], j)
            cat(text, file=log_file0, append = TRUE)
            if (ignore == FALSE){
              cat(text)
              if(nP$ALL.removed){
                cat('Do you want to abort the process?\nPlease enter a number from the selection\n1 = Yes\n2 = No\n3 = Ignore further warnings\n')
                interact <- readLines('stdin', n=1)
                if(interact == 1){
                  stop('Pipeline stopped')
                } else if(interact == 2){
                  cat('Found primers/MIDs were written to config file')
                } else {
                  ignore <- TRUE
                }
              }
            }           
          }           
        }
    } 
  }
  
}
 

##### Step 1 #####
read.Assembler <- function(input_folder, threshold, minoverlap, minqual, minlength, maxlength, primer.offset, pre.filter, pre.filter.mq){

  ## uses pandaseq to assemble paired End reads  
  sys = system("pandaseq -v 2>&1 ", intern = T)
  cat(sys[1], "\n\n", file = log_file1 , append = T)

  
  ############# defines used primer sequence #########################################################
  ## poly-N, MID, specific forward primer
  primer  <- read.csv(nP$primerTable, header=T, stringsAsFactors = F)
  comb.f = c() ## forward primer
  comb.r = c() ## reverse primer
  ## pandaseq just takes the length of the primers but not the sequence, in case primer are different
  ## and can' be expressed by common 1L code
  if(primer.offset == TRUE){
    # cat(paste(Sys.time(), "| function read.Assembler with primer.offset=TRUE \n", sep=""), file = nP$log_file , append = T)
    if (nP$MID.removed){ ## only forward primer
      # cat(paste(Sys.time(), "| function read.Assembler with 'MID.removed' \n", sep=""), file = nP$log_file , append = T)
      for(i in 1:length(primer[,1])){
        tmp = primer[i,"specific_forward_primer"]
        tmp = unlist(strsplit(tmp, split = ""))
        comb.f <- c(comb.f, length(tmp))
        tmp = paste0(primer[i, "poly_N_rev"], primer[i, "specific_reverse_primer"], collapse = "")
        tmp = unlist(strsplit(tmp, split = ""))
        comb.r <- c(comb.r, length(tmp))
      }
    }else{ ### takes primer sequence
      # cat(paste(Sys.time(), "| function read.Assembler without 'MID.removed' \n", sep=""), file = nP$log_file , append = T)
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
    # cat(paste(Sys.time(), "| function read.Assembler with primer.offset=FALSE \n", sep=""), file = nP$log_file , append = T)
    if (nP$MID.removed){
      # cat(paste(Sys.time(), "| function read.Assembler with 'MID.removed' \n", sep=""), file = nP$log_file , append = T)
      for(i in 1:length(primer[,1])){
        comb.f <- c(comb.f,  as.character(primer[i,"specific_forward_primer"]))
        comb.r <- c(comb.r, as.character(primer[i, "specific_reverse_primer"]))
      }
    }else{
      # cat(paste(Sys.time(), "| function read.Assembler without 'MID.removed' \n", sep=""), file = nP$log_file , append = T)
      for(i in 1:length(primer[,1])){
        comb.f <- c(comb.f, paste( primer[i,"MID"], primer[i,"specific_forward_primer"], sep=""))
        comb.r <- c(comb.r, as.character(primer[i, "specific_reverse_primer"]))
      }
    }
    primer <- cbind(primer, f_primer = comb.f, r_primer = comb.r)
  }
  ##################################################################################################
  ########################  iterates over all files    ############################################
  files_forward = list.files(input_folder, full.names=TRUE, pattern="\\R1.fastq$")
  files_reverse = list.files(input_folder, full.names=TRUE, pattern="\\R2.fastq$")
  
  prinseq <- "./bin/prinseq-lite/prinseq-lite.pl"
  cat(system(paste(prinseq, " --version", sep = ""), intern = T), "\n\n", file = log_file1, append =T)
  
  for(i in (1:length(files_forward))){
    time = date()
    
    file.f <- files_forward[[i]]
    file.name <- unlist(strsplit(file.f, split = "/"))
    file.name <- file.name[length(file.name)]
    file.name <- unlist(strsplit(file.name, split = "_R1.fastq" ))
    
    file.r <- files_reverse[[i]]
    file.name.r <- unlist(strsplit(file.r, split = "/"))
    file.name.r <- file.name.r[length(file.name.r)]
    file.name.r <- unlist(strsplit(file.name.r, split = "_R2.fastq" ))
    
    folder.sys.call = system(paste("mkdir -p ",nP$resultsFolder, "/", file.name , sep=""), intern = TRUE)
    localFolder <- paste(nP$resultsFolder, "/", file.name , sep="")
    log_folder <- paste(input_folder, "/logging/assembling", sep="")
    
    ### remove reads with an average quality below treshold
    if(pre.filter == TRUE){
      cat(paste(Sys.time(), "| start quality filtering \n", sep=""), file = log_file1 , append = T)
      out.good <- paste(localFolder, "/", file.name, sep ="")
      out.bad <- paste(localFolder, "/", file.name, "_bad", sep = "")
      filter.cmd = paste(prinseq, " -verbose -fastq ", file.f, " -fastq2 " , file.r, " -ns_max_n 0 -min_qual_mean ", pre.filter.mq, " -out_good ", 
                         out.good ," -out_bad ", out.bad, " 2>&1 " , sep = "")
      sys = system(filter.cmd, intern = T) 
      sys = sys[- c(1:6)]
      cat(paste(sys, sep = " " , collapse = "\n"), "\n\n", file = log_file1 , append = T)
      
      file.f <- paste(out.good, "_1.fastq", sep = "")
      file.r <- paste(out.good, "_2.fastq", sep = "")
    }
    
    cat(paste(Sys.time(), "| Start assembling ", file.name, "\n", sep = ""), file = log_file1 , append = T)
    cat("read.Assembler: \n")
    if(file.name == file.name.r){
      this.primer <- primer[primer$Probe_FWD == file.name,]
      this.primer.f <- as.character(this.primer$f_primer)
      this.primer.r <- as.character(this.primer$r_primer)
      
      logging.file <- paste(log_folder ,"/pandaseq_log.txt", sep = "")
      panda.out <- paste(localFolder, "/", file.name, "_assembled.fastq", sep ="")
      
      cat("Checking for existence of some files in ", localFolder, ": \n")
      cat("Forward: ", file.exists(file.f), " ",file.f, "\n")
      cat("Reverse: ", file.exists(file.r), " ", file.r, "\n") 
      
      ## call PandaSeq
      if(nP$ALL.removed){
        cmd = paste("pandaseq -f ", file.f, " -r ", file.r, " -B -a -F -g ", logging.file, " -w ", panda.out, " -N  -t ", threshold, " -o ", minoverlap, " -l ", minlength,
                    " -L ", maxlength, " -C min_phred:", minqual , sep= "" )
      }else{
        cmd = paste("pandaseq -f ", file.f, " -r ", file.r, " -B -a -F -g ", logging.file, " -w ", panda.out, " -N -p \"",
                    this.primer.f, "\" -q \"", this.primer.r, "\" -t ", threshold, " -o ", minoverlap, " -l ", minlength,
                    " -L ", maxlength, " -C min_phred:", minqual , sep= "" )
      }
      
      sys = system(cmd, intern = T) 
      cat(paste(Sys.time(), "| Statement executed: \n", cmd, "\n\n", sep=""), file = log_file1, append = T)
      
      cat("panda logging: ", file.exists(logging.file)," ", logging.file, "\n")
      cat("panda assembled: ", file.exists(panda.out), " ", panda.out, "\n\n")
      
      ### read fastq file, rename sequences and save as fasta
      if(file.info(panda.out)$size > 0){
        assembled.fastq <- readFastq(panda.out)
        assembled.seq <- sread(assembled.fastq)
        
        nr <- length(assembled.fastq)
        cat(paste(Sys.time(), "| assembled ", nr, " reads \n \n", sep=""), file = log_file1 , append = T)
        new.names <- sapply(X = 1: nr, FUN = function(x, file.name){
          name <- paste(file.name, "_", x, sep= "")
          return(name)
        }, file.name = file.name)
        
        names(assembled.seq) <- new.names
        writeXStringSet(assembled.seq, filepath = paste(localFolder, "/", file.name , ".fasta", sep=""))
        file.create(paste0(log_folder1, '/', 'assembling_success'))

      }else{ ## no reads assembled
        main.folder <- strsplit(nP$resultsFolder, split = "/results")[1] 
        ## make folder for sample with no assembled reads
        bad.assembled.folder <- paste(main.folder, "/not.assembled", sep ="")
        if(! file.exists(bad.assembled.folder)){
          cmd = paste("mkdir -p ", bad.assembled.folder, sep ="" )
          sys = system(cmd, intern = T)   
        }        
        ## move folder
        cmd = paste("mv", localFolder,  bad.assembled.folder, sep =" ")
        sys = system(cmd, intern = T)
        message = "no reads assembled \n\n" 
        cat(message, file = log_file1 , append = T)
        
      }
      
    }else{
      stop("forward and reverse file don't match!")
    }     
  } 
}


##### Step 2 #####
cdhit.clustering100 <- function(input, output, cores, th){
  # cat("Start of function cdhit.clustering100 \n")
  cat(paste("Variables filled with: \n input =", input, "\n output = ", output, "\n cores = ", cores, "\n th = ",th, "\n"))
  #### functions ####
  getlabels <- function(name, frequencies){
    labels <- sapply(X = 1:length(frequencies), FUN = function(x){
      label <- paste(name,"_", x,";size=", frequencies[x], ";", sep = "")
      return(label)
    })
  }
  if(file.exists(input)){
    ### run cd-hit ####
    representative.path = sub(".fasta", "_cdhit", input)
    cluster.path=sub(".fasta", "_cdhit.clstr", input)
    cdhit.path = './bin/cd-hit-v4.6.5-2016-0304'
    command <- paste(cdhit.path, "/cd-hit-est -i ./",input, " -o ./", representative.path, " -c ", th, " -T " , cores, sep = "")
    syst = system(command, intern = T)
    cat("command: ", command, "\n\n")
    
    #### get cluster sizes ###########
    #read cluster file
    representatives <- readDNAStringSet(representative.path)
    #cluster <- scan(cluster.path, what= "character", )
    cluster <- readLines(cluster.path)
    ## number of reads
    cluster.nr <- length(representatives)
    read.nr <- length(cluster)- cluster.nr
    ## start position new cluster in cluster list
    match <- grep(">Cluster", cluster)
    ## get size of each cluster
    match.shift <- as.integer(match + 1)
    match <- c(match[-1],length(cluster)+1 )
    cluster.length <- match - match.shift
    cat(paste(Sys.time(), "| ", input, ":\n number of reads sorted in clusters => ",sum(cluster.length), "\n\n", sep=""), file = log_file2 , append = T)
    
    ## rename representative cluster, add count information
    cluster.names <- getlabels("OTU", cluster.length)
    names(representatives) <- cluster.names
    
    ## arrange cluster by decreasing size  
    desc.order <- order(cluster.length, decreasing = TRUE)
    representatives <- representatives[desc.order]
    
    #### save fasta with size information
    writeFasta(representatives, output)
    cat("File has been written: ", file.exists(output), " --> ",  output, "\n")
    # cat("End of function cdhit.clustering100 \n\n")
    rm(match, match.shift, cluster, representatives)
    gc()
    return(syst)
  }else{
    #### TODO check if stop
    cat("File ", input, "does not exist. \n\n")
  }
}


##### Step 3 optional #####
getPositiveChimeras <- function(uchimeout.file, input.file, nonchimera.file){
  #cat("Function getPositiveChimeras \n")
  cat("uchimeout.file = ", uchimeout.file, "\n")
  cat("input.file = ", input.file, "\n")
  cat("nonchimera.file = ", nonchimera.file, "\n\n")
  ## read input fasta
  fasta <- readFasta(input.file)
  fasta.ids <- as.character(id(fasta))
  ## reads uchime results
  uchime.table <- read.table(uchimeout.file,stringsAsFactors=FALSE)
  positives <- which(uchime.table$V18 == "Y")
  if(length(positives) > 0){ ## chimeras found
    uchime.table.pos <- uchime.table[positives, ]
    chimeras.names <- uchime.table.pos$V2
    ### find chimera positions in input fasta
    chim <-which( is.element(fasta.ids, chimeras.names))
    ## remove chimeras
    fasta.wo.chimeras <- fasta[- chim]
    ## write new fasta file wo chimeric sequences
    writeFasta(fasta.wo.chimeras, file = nonchimera.file)
    cat("chimera length >0; written nonchimera.file: ", file.exists(nonchimera.file), "\n ", nonchimera.file, "\n\n")
  }else{## no chimeras found
    writeFasta(fasta, file = nonchimera.file)
    cat("written nonchimera.file: ", file.exists(nonchimera.file), "\n ", nonchimera.file, "\n\n")
  }
  
  return(NULL)
}

######### function definition: iterate over all cluster and counts read  ############################
count.cluster.size <- function(cluster.dir, saving.dir, sample.name){
  
  cluster <- read.fasta(cluster.dir) ## fasta files of cluster
  centroid <- cluster[[1]] # list of one SeqFastadna object
  
  names.members <- names(cluster)
  size.members <-  sapply(X = names.members, FUN= function(x){
    tmp = unlist(strsplit(x, split=";size="))[2]
    tmp = unlist(strsplit(tmp, split=";"))[1]
    return(as.integer(tmp))
  })
  
  cluster.size <- sum(size.members)
  new.name <- unlist(strsplit(getAnnot(centroid), split=";size="))
  cluster.name <- new.name[1]
  new.name[2]<- paste(cluster.size, ";", sep = "")
  new.name <- paste(new.name[1],";size=", new.name[2], sep = "")
  centroid.seq = unlist(getSequence(centroid, as.string=T))
  save.name <-  paste(saving.dir, "/", sample.name, "_","centroids.abundance.fasta", sep = "")
  
  sequence=unlist(getSequence(centroid.seq, as.string=T))
  upper.seq <- sapply(sequence, toupper)
  
  write.fasta(upper.seq, names = new.name, file.out =  save.name, open = "a")
  
  #save.name1 <-  paste(saving.dir, "/", sample.name, ".clustered", similarity, ".fasta", sep = "")
  #cmd <- paste("./bin/usearch7 -sortbysize ", save.name, "-output ",  save.name1, " -minsize 1")
  #sys = system(cmd, intern = T)
  return(cluster.size)
}
### end function count.cluster.size ##########################################################


##### Step 4 #####

# Important:
# + the user must make sure that the correct paths for the input is provided.
# + in case the package "seqinr" is not pre-installed, an internet connection 
#   is needed. 
# + in some cases it was necessary to manually install the package Rcpp, in order for the dependent
#   package data.table to be installed.
#
#
# Author: 
# Simo Kitanovski, simo.kitanovski@stud.uni-due.de, 02.05.2014

fastaPoolReader <- function(pathToFastaFiles) {
  cat("\n Start of function fastaPoolReader \n")
  
  # function fastaPoolReader - it takes a path as an input at which all the FASTA files
  # are located. It then reads them one by one and creates a single csv file with the data 
  # extracted from each file. An output of 2 FASTA files would have the following column 
  # structure: sequences | fileA | fileB whereby the colums fileA and fileB would contain 
  # the size parameter for each sequence. This file is created at dir output/mergedFasta.csv,
  # which is created if it doesn't exist.
  #  
  #Check if seqinr can be loaded, else try to install it. 
  #Then try loading it again.
  # require(seqinr),is already loaded in the pipeline 
  
  #input csv pool
  fastaFiles <- list.files(path = pathToFastaFiles, pattern = "*.fasta")
  nrFastaFiles <- length(fastaFiles)
  
  sample.names <- sapply(X = fastaFiles, FUN= function(x){
    tmp <- unlist(strsplit(x, split = ".", fixe = T))[1]
    return(tmp)
  })
  
  #name/dir of the output csv file. If the output directory doesnt exist, create it.
  outDir <- paste(strsplit(pathToFastaFiles, split = "/finalData")[1], "/tables/mergedFasta.csv", sep ="")
  cat("checking the variables: \n fastaFiles: ", fastaFiles, "\n nrFastaFiles: ", nrFastaFiles, "\n outDir: ", outDir, "\n\n")
  #if the file mergedFasta.csv exists, delete it.
  if(file.exists(outDir)) {
    file.remove(outDir)
  }
  
  #loop through the FASTA files, read them one by one as csv files.
  for(file in seq(from = 1, to = nrFastaFiles, by = 1)) {
    
    #get the name of the FASTA file (the entire directory included).
    currentFile <- paste(pathToFastaFiles, "/", fastaFiles[file], sep = "")
    
    #read the fasta file, print a message which file is currently being read.
    cat(paste("merging FASTA file:", fastaFiles[file] , format(Sys.time(), "%R"), "\n", sep = " "))
    cat(paste("merging FASTA file:", fastaFiles[file] , format(Sys.time(), "%R"), "\n", sep = " "), file = log_file4 , append = TRUE)
    fastaFile <- read.fasta(file = currentFile, as.string = T)
    
    #obtain the sequences only from the FASTA file.
    sequences <- getSequence(object = fastaFile, as.string = T)
    
    #unlist the sequences, creating a data frame with 1 column where they are stored.
    dfSequences <- data.frame(sequences = unlist(sequences))
    
    #obtain the annotations only from the FASTA file.
    annotations <- getAnnot(object = fastaFile, as.string = T)
    
    #the next three commands are used to obtain the 'size' parameter for each sequence.
    #The number variable is the final result which contains all the 'sizes' in the same
    #order as the sequences are stored in dfSequences.
    size <- unlist(annotations)
    size <- sapply(strsplit(size, ";"), as.character)
    #this is now a column of sizes.
    number <- sapply(strsplit(size[2, ], "="), as.character)
    rm(size)
    #loop through the FASTA file names. For each file name do one of the following two options:
    # + IF the file name in this inner loop is not equal to the file name of the outer loop (file.in != file)
    #   => append a column to the dfSequences data.frame with all o's.
    # + ELSE (file.in == file)  
    #   => append a column with the actual 'sizes' stored previously in the variable 'number' to the dfSequences 
    #   data.frame
    
    for(file.in in seq(from = 1, to = nrFastaFiles, by = 1)) {
      if(file.in == file) {
        dfSequences <- cbind(dfSequences, number[2, ])
      }
      else {
        dfSequences <- cbind(dfSequences, 0)
      }
    } 
    #export the final data.frame as a csv file to the output dir.
    colnames(dfSequences) <- c("sequence", sample.names)
    ## appends new table to previous table
    write.table(x = dfSequences, file = outDir, append = T, quote = F, col.names = F,row.names = F, sep=";")
    cat("New table has been written: ", file.exists(outDir), outDir, "\n")
  }  
  cat(paste(Sys.time(), "| All files combined in one table: ", outDir, "\n", sep = ""), file = log_file4 , append = TRUE)
  cat("End of function fastaPoolReader \n\n")
}


hashBuilder <- function(pathToFiles, table.name) {
  cat("Start of function hashBuilder \n")
  # function hashBuilder - it takes the output of the output created by the previous function 
  # fileReader. For any two duplicated sequence it adds up their corresponding values from each 
  # other column. Finally the hash table is exported as a csv at dir output/hashtable.csv, which
  # is created if it doesn't exist.
  #
  #Check if data.table can be loaded, else try to install it. 
  #Then try loading it again.
  if(require("data.table") == F) {
    install.packages("data.table") 
    require("data.table")
  }
  
  #name/dir of the output csv file. If the output directory doesnt exist, create it.
  if(nP$filtering){
    outDir <-  paste(pathToFiles,"/", table.name, "_unfilteredTable.csv", sep="")
  }else{
    outDir <-  paste(pathToFiles,"/", table.name, "_Table.csv", sep="")
  }
  cat("This outDir will be used: ", outDir, "\n")
  pathToDir <- paste(strsplit(pathToFiles, split = "/finalData")[1], "/tables/mergedFasta.csv", sep ="")
  cat("And this is the pathtoDir: ", pathToDir, "\n")
  #list all the files in the input directory (pathToFiles)
  fastaFiles <- list.files(path = pathToFiles, pattern = "*.fasta")
  cat("These fastaFiles will be used: ", fastaFiles, "\n")
  sample.names <- sapply(X = fastaFiles, FUN= function(x){
    tmp <- unlist(strsplit(x, split = ".", fixe = T))[1]
    return(tmp)
  })
  
  #do fast csv reading (provided by data.table package). It reads as input the csv file 
  #produced by the fastaPoolReader function. Print a message with the system time when 
  #the process starts/ends.
  dt <- fread(input = pathToDir, sep = ";", header = F, stringsAsFactors = F)
  #sets the names of the data table. Identical as in the colnames function.
  setnames(dt, old = colnames(dt), new = c("sequences", sample.names))
  #setattr(dt, name = colnames(dt), value = c("sequences", sample.names))
  
  #aggregation = wherever the sequences are duplicate, sum up the frequencies included 
  #in the other columns.
  agg <- dt[, lapply(.SD, sum), by = sequences]
  
  tables.info <- tables()
  text <- paste("initial table has ", tables.info[,NROW][2], " rows and occupies ", tables.info[,MB][2], " MB memory. \n", sep = "")
  cat(text, file = log_file4 , append = TRUE)
  text <- paste("aggregated table has ", tables.info[,NROW][1], " rows and occupies ", tables.info[,MB][1], " MB memory. \n", sep = "")
  cat(text, file = log_file4 , append = TRUE)
  
  #write final output file as a csv
  write.table(x = agg, file = outDir, row.names = FALSE)  
  cat("File has been written: ", file.exists(outDir), outDir, "\n")
  cat("typeof(agg): ", typeof(agg), "\n")
  cat("End of function hashBuilder \n\n")
  return(as.data.frame(agg))
}



##################################### make table function definition   ###############
#The make.table function is replaced by two functions fastaPoolReader and hashBuilder
#included in the script Hasher.R

make.table <- function(data.folder, table.name){
  cat("\n\n", paste(Sys.time(), "| Start making count table: \n", sep=""), file = log_file4 , append = TRUE)
  cat("data.folder: ", data.folder, "\n")
  cat("table.name: ", table.name, "\n")
  # source("./bin/Hasher.R")

  table <- fastaPoolReader(pathToFastaFiles = data.folder)

  table = hashBuilder(pathToFiles = data.folder, table.name)
  return(table)
}


################################################################################################################################################
################################################### singelton.filter ###########################################################################
################################################################################################################################################
# singelton -> only one occurance over all samples
# cutoff min number of reads
# min number of samples that should have the sequence
singelton.filter <- function(table, final.data.folder, cutoff, occurance){
  cat(paste( "\n \n", Sys.time(), "| start singelton filtering. \n", sep=""), file = log_file4 , append = TRUE )
  
  filtered = table
  sample.count = length(names(filtered)) -1
  
  cat(paste( "\n \n These parameters are provided for singleton filtering: \n",sep = ""), file = log_file4 , append = TRUE )
  if(!missing(table)){
    cat("The variable 'table' is filled \n")
  }else{cat("'table' is not filled")}
  if(!missing(final.data.folder)){
    cat("The variable 'final.data.folder' is filled \n")
  }else{cat("'final.data.folder' is not filled")}
  if(!missing(cutoff)){
    cat("The variable 'cutoff' is filled \n")
  }else{cat("'cutoff' is not filled")}
  if(!missing(occurance)){
    cat("The variable 'occurance' is filled \n")
  }else{cat("'occurance' is not filled")}
  if(!missing(filtered)){
    cat("The variable 'filtered' is filled \n")
  }else{cat("'filtered' is not filled")}
  if(!missing(sample.count)){
    cat("The variable 'sample.count' is filled \n")
  }else{cat("'sample.count' is not filled")}
  
  
  ## sets all frequencies below the cutoff to zero
  tmp <- foreach(i = 1:length(table$sequences), .combine="c") %do% {
    occ <- sum(filtered[i, -1 ] > 0)
    cut <- sum(filtered[i, -1 ])
    if( (occ >= occurance) && (cut >= cutoff)){
      value=T
    }else{
      value=F
    }
    value
  }
  
  obsolete_sequences = !tmp
  ## removes all sequences that have no occurences after filtering
  obsolete_sequences.table <- table[obsolete_sequences,]
  outDir.f <- paste(final.data.folder,"/",  nP$input_folder, "_filtering_failed.csv", sep="")
  write.table(x = obsolete_sequences.table, file = outDir.f, row.names = FALSE)
  rm(obsolete_sequences.table)
  
  if(length(obsolete_sequences) > 0){
    filtered <- filtered[tmp, ]
  }
  
  cat(paste("\n kept ", length(filtered[,1]) ," sequence \n \n", sep = ""), file = log_file4 , append = TRUE )
  
  ## saves filtered data table
  outDir <-  paste(final.data.folder,"/",  nP$input_folder, "_Table.csv", sep="")
  write.table(x = filtered, file = outDir, row.names = FALSE)
  cat("filtered data has been written: ", file.exists(outDir), outDir, "\n")
  cat(paste( "\n \n", Sys.time(), "| End of singelton filtering. \n \n", sep = ""), file = log_file4 , append = TRUE )
}



#################################################################################################################################################
########################################################  mixed.filtering  ######################################################################
#################################################################################################################################################

mixed.filtering <- function(table, final.data.folder, cutoff, p.correction, saving.format, nP){
  cat(paste( "\n \n", Sys.time()," start mixed filtering. \n \n", sep = ""), file = log_file4 , append = TRUE )
  
  figure.folder = paste(final.data.folder, "/figures", sep="")
  sys = system(paste("mkdir ", figure.folder, sep=""), intern = TRUE)
  
  
  if(nP$pairedEnd == FALSE){
    ext <- nP$name_extension
    sample.names <- names(table)
    sample.names <- gsub(ext, "" , sample.names)
    names(table) <- sample.names
  }
  
  sample.names <- names(table)
  sample.names <- sample.names[-1] # +1 offset from table
  sample.count =  length(sample.names)
  patternA = "*.A$"
  matchA = grep(patternA, sample.names)
  patternB = "*.B$"
  matchB = grep(patternB, sample.names)
  
  individual.samples = c()
  
  ## looks for samples, that are neither A nor B
  if((length(matchB) + length(matchA)) != sample.count){
    AB_samples = c(matchA, matchB)
    individual.samples <- c(1:sample.count)
    individual.samples <- setdiff(individual.samples, AB_samples)
  }
  
  ## looks if every A sample has a B sample
  ## A-B paires are stored in AB_samples
  ## A samples that have no B sample are appenden to individual samples
  AB_samples = c()
  res <- foreach(i = 1:length(matchA)) %do% {
    j = matchA[i]
    sample <- unlist(strsplit(sample.names[j], split = "A$"))
    tmp <- grep(paste("^", sample, "B$", sep = ""), sample.names[matchB])
    if(length(tmp) == 1){
      hit = matchB[tmp]
      AB_samples <- c(AB_samples, j, hit)
    }else{individual.samples = c(individual.samples, j)
    }  
  }
  
  ## looks for B samples, that have no corresponding A sample
  individualBs <- setdiff(matchB, AB_samples)
  individual.samples <- c(individual.samples, individualBs)## numbering without sequence column, offset +1
  
  filtered <- table
  
  ## ampliconduo before filtering
  if(length(AB_samples) >= 2){
    require(AmpliconDuo)
    intable <- table[, -1]
    intable <- intable[, AB_samples]
    names <- names(intable)
    names <- names[seq(1, length(names), by = 2)]
    names <- sapply(X = names, FUN = function(x){
      unlist(strsplit(x, split = "+.A$"))
    })
    ampliconduos <- ampliconduo(intable, sample.names = names, correction = p.correction) 
    rm(intable)
    
    if(nP$plot.ampliconduo){
      nrow = as.integer(sqrt(length(names)))
      tmp <-plotAmpliconduo.set(ampliconduos, nrow = nrow, save = TRUE, path = figure.folder, format = saving.format)
    }
    
    ## discordance Delta befor filtering
    tmp <- discordance.delta(ampliconduos, corrected = TRUE, printToTex= TRUE, directory = figure.folder)
    ## filtering #############
    
    
    ## apply ampliconduo filtering only to AB_samples
    tmp <- foreach(i = 1:(length(AB_samples)/2)) %do% {
      index_a <- AB_samples[(2*i -1)] + 1
      index_b <- AB_samples[(2*i)] +1
      bad.amps <- which(filtered[,index_a] < 1| filtered[,index_b] < 1)
      if(length(bad.amps) > 0){
        filtered[bad.amps, index_a] <- 0
        filtered[bad.amps, index_b] <- 0
      }     
    }
  }
  
  ### apply cutoff filter to individual samples
  if(length(individual.samples) > 0){
    tmp <- foreach( i = 1:length(individual.samples)) %do% {
      index <- individual.samples[i] +1
      bad.amps <- which(filtered[, index] < cutoff)
      if(length(bad.amps) > 0){
        filtered[bad.amps, index] <- 0
      }
    } 
  }
  
  
  ### remove all rows with only zeros ###########
  if(sample.count > 1){
    obsolete_sequences <- which(rowSums(filtered[, 2:(1+sample.count)]) == 0)
  }else{
    obsolete_sequences <- which(filtered[, 2] == 0)
  }
  if(length(obsolete_sequences) > 0){
    obsolete_sequences.table <- table[obsolete_sequences,]
    outDir.f <- paste(final.data.folder,"/",  nP$input_folder, "_filtering_failed.csv", sep="")
    write.table(x = obsolete_sequences.table, file = outDir.f, row.names = FALSE)
    rm(obsolete_sequences.table)
    
    filtered <- filtered[- obsolete_sequences, ]
  }
  
  
  ## ampliconduo with filtered samples
  if(length(AB_samples) >= 2){
    intable <- filtered[, -1 ]
    ## find samples with no reads left
    read.counts <- colSums(intable)
    zero.reads.count <- which(read.counts == 0)
    if(length(zero.reads.count) > 0){
      zero.reads.count.names <- names(intable[zero.reads.count])
      AB_samples.sub.zero.read.samples <- setdiff(AB_samples, zero.reads.count)
      zero.reads.count.names = paste0(zero.reads.count.names, collapse = "\n")
      cat(paste( "\nSamples with 0 reads after filtering removed:\n", zero.reads.count.names ,"  \n \n", sep = ""), 
          file = log_file4 , append = TRUE )
    }else{
      AB_samples.sub.zero.read.samples <- AB_samples
    }
    
    intable <- intable[, AB_samples.sub.zero.read.samples]
    names.intable <- names(intable)
    names.intable <- names[seq(1, length(names.intable), by = 2)]
    names.intable <- sapply(X = names.intable, FUN = function(x){
      unlist(strsplit(x, split = "+.A$"))
    })
    ampliconduos.f <- ampliconduo(intable, sample.names = names.intable, correction = p.correction) 
    rm(intable) 
    
    ### plot ampliconduo
    if(nP$plot.ampliconduo){
      file.name <- paste("ampliconduo_filtered_", Sys.Date(), sep = "")
      tmp <-plotAmpliconduo.set(ampliconduos.f, nrow = nrow, save = TRUE, path = figure.folder, format = saving.format
                                , file.name = file.name)
    }
    
    ### discordance delta
    file.name <- paste("discordanceDelta_filtered_", Sys.Date(), sep = "")
    tmp <- discordance.delta(ampliconduos.f, corrected = TRUE, printToTex= TRUE, directory = figure.folder, file.name = file.name) 
    
    save(ampliconduos, ampliconduos.f, file = paste(figure.folder, "/AmpliconDuo.RData", sep = ""))
  }
  if(length(zero.reads.count) > 0){
    filtered <- filtered[, - (zero.reads.count + 1)]
  }
  
  cat(paste( "\nkept ", length(filtered[,1]) ," sequence \n \n", sep = ""), file = log_file4 , append = TRUE )
  
  outDir <-  paste(final.data.folder,"/",  nP$input_folder, "_Table.csv", sep="")
  write.table(x = filtered, file = outDir, row.names = FALSE)
}

swarm <- function(final.data.folder, cores){
  home.dir <- getwd()
  setwd(final.data.folder)
  swarm <- paste0('swarm -t ',cores, ' -f -z -w merged.representatives.fasta <  merged.fasta > merged.swarms')
  system(swarm)
  setwd(home.dir)
}
