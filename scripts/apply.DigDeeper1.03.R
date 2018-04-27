#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| apply.DigDeeper1.03.R                               |\n")
#cat("|                                                     |\n")
#cat("| (c) by Anja Lange, 2014                             |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("| anja.lange@uni-due.de                               |\n")
#cat("| implements functions by Dominik Heider              |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")

# The script was kept mostly as-is, I just added the function call at the end of it -Marius.
# have to adjust a few more things, get rid of the global env crap and finish the function call at the end of the script
library('ShortRead')
with_primertable = TRUE
primertable = read.csv('p_table_with_primer', header = TRUE, stringsAsFactors= FALSE, sep = ",")

apply.DigDeeper <- function(file, sample.name, primer, max.length){
  ## directory + file name file
  ## sample.name is the name with the extension like _R1
  ## make results folder for sample
  folder.sys.call = system(paste("mkdir ",.GlobalEnv[["resultsFolder"]], "/", sample.name , sep=""), intern = TRUE)
  localFolder <- paste(.GlobalEnv[["resultsFolder"]], "/", sample.name , sep="")
  
  ## read data
  sample <- readFastq(file) ## ShortRead
  local.log <- paste(Sys.time(), "\n", file, "\n", "loaded ", length(sample) , " sequences \n", sep = "")
  
  #################################################################################
  #### Poly N trimming
  trim_results <- trim(sample, sample.name, local.log, localFolder)
  local.log <- trim_results[[1]]
  trim_results <- trim_results[[2]]
  
  ################################################################################
  ### remove short sequences
  
  if (min_length_of_sequence){
    first_removed <- remove_sequences_below_lengthcutoff(trim_results, sample.name, primer, .GlobalEnv[["with_primertable"]],
                                                         .GlobalEnv[["min_length_of_sequence"]], local.log, localFolder)
    local.log <- first_removed[[1]]
    first_removed <- first_removed[[2]]
    
  }else{
    first_removed <- trim_results
  }
  
  ####################################################################################
  ##### check quality
  # Pakete der Größe 250.000 bauen, weil Illumina sonst den Speicher voll macht 
  # ceiling takes a single numeric argument x and returns a numeric vector containing
  # the smallest integers not less than the corresponding elements of x.
  num.packages = ceiling(length(first_removed) / 250000)
  for(i in 1:num.packages ){
    start = 250000 * i - 250000 + 1
    stop = 250000 * i
    if(stop>length(first_removed)){
      stop = length(first_removed)
    }
    cat(paste("Package number: ", i, "/", num.packages, "\n", sep=""))
    local.log <- paste(local.log, "Package number: ", i, "/", num.packages, "\n", sep="" )
    quality_results <- qualy(first_removed[start:stop], .GlobalEnv[["cut_off_base"]], .GlobalEnv[["cut_off_seq"]],
                             FALSE, sample.name, local.log, localFolder)
    local.log <- quality_results[[1]]
    print(start)
    print(stop)
  }
  
  quality_results <- readDNAStringSet(paste(localFolder, "/", sample.name , "-checked.fasta", sep=""))
  cat("In total ", length(quality_results), " sequences passed Quality control \n" , sep = "")
  local.log <- paste(local.log, "In total ", length(quality_results), " sequences passed Quality control \n" , sep = "")
  ####################################################################################
  ###### remove forward primer
  
  #quality_results are now in fasta format! Not in fastq
  if (with_primertable){  
    if(.GlobalEnv[["ALL.removed"]]){ ## no primer sequence that needs to be removed, just rename quality results and save it
      nr <- length(quality_results)
      new.names <- sapply(X = 1: nr, FUN = function(x, samplename){
        name <- paste(samplename, "_", x, sep= "")
        return(name)
      }, samplename= sample.name)
      
      names(quality_results) <- new.names
      writeXStringSet(quality_results, filepath = paste(localFolder, "/", sample.name , ".fasta", sep=""))
    }else{
      ## removed contains all the sequences that were identified using the MID and primer, and these sequence parts are removed
      ## the sequence names correspond now to the sample
      ## removed <- removeForwardPrimer(quality_results, primerTable, verbose, filename, multicores)
      all_data <- foreach (p = 1: length(quality_results), .combine="c") %dopar% {
        one_seq <- primer_cutter(quality_results[p], primer, max.length)
      } 
      
      removed <- all_data ### DNAStringSet with removed primers    
      writeXStringSet(removed, filepath = paste(localFolder, "/", sample.name , "-prim_removed.fasta", sep=""))
      
      number_mismatched <- length(quality_results) - length(removed)
      local.log <- paste(local.log, "forward primer removed: ", date(), "\n", number_mismatched, 
                         " sequences did not match the primer \n" , sep = "")
      
      cat("Forward primer removed.\n \n")
      
      ## removes short sequences that are now below some treshhold
      ## obsolete?
      result <- (removed[width(removed)>=min_length_of_sequence])
      writeXStringSet(removed, filepath = paste(localFolder, "/", sample.name , "-wo_ShortSeq2.fasta", sep=""))
      
      ### rename the sequences sample name + number
      nr <- length(result)
      new.names <- sapply(X = 1: nr, FUN = function(x, samplename){
        name <- paste(samplename, "_", x, sep= "")
        return(name)
      }, samplename= sample.name)
      
      names(result) <- new.names
      writeXStringSet(result, filepath = paste(localFolder, "/", sample.name , ".fasta", sep=""))
    }    
  } 
  
  ######################################################################################  
  #### append local.log to the logfile
  cat(local.log, file = .GlobalEnv[["logfile"]] , append = TRUE)
  cat("\n \n", file = .GlobalEnv[["logfile"]] , append = TRUE)
}

##########
#####################################################################################################################

#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| Poly-N trimmer v0.1                                 |\n")
#cat("|                                                     |\n")
#cat("| (c) by Dominik Heider, 2012                         |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("| dominik.heider@uni-due.de                           |\n")
#cat("| adapted by Anja Lange, 2014                         |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")

trim <- 
  function(data, samplename, local.logfile, localFolder){
    
    cat(paste(samplename, ": loading Fastq file...\n"))
    reads <- data
    logfile <- local.logfile
    logfile <- paste(logfile, "Before Trimming: ", length(reads), " sequences \n")
    #######################################
    ## Inject Ns at low quality positions ## part was removed!! filters low quality sequences later
    ########################################
    seqs <- sread(reads) #define ‘accessors’ (to get and set values) for objects in the ShortRead package  get sequence list
    qual <- PhredQuality(quality(quality(reads))) #first quality() converts input to FastqQuality, second quality()
    ### to a BStringSet and finally to a PhredQuality instance
    injectedseqs <- seqs
    
    cat("removing poly-N tails...\n")
    
    ####################################
    ## Get coordinates of polyN tails ##
    ####################################
    adapter <- paste(rep("N", max(width(injectedseqs))), sep="", collapse="") #String of Ns as pattern to identify polyN tails
    mismatchVector <- c(rep(0,width(adapter))) # vector of 0s ,allow no mismatches at each adapter offset
    # The trimLRPatterns function trims right flanking patterns from sequences. gives the coordinates of the remaining,
    # trimmed function back
    trimCoords <- trimLRPatterns(Rpattern=adapter, subject=injectedseqs, max.Rmismatch=mismatchVector, ranges=T)
    # Trim sequences looking for a right end pattern (polyN in this case)
    # Gets IRanges object with trimmed coordinates
    
    cat("saving trimmed file...\n")
    ##########################################################################
    ## Apply trimming coordinates from injected reads to non-injected reads ##
    ##########################################################################
    seqs <- DNAStringSet(seqs, start=start(trimCoords), end=end(trimCoords))
    qual <- BStringSet(qual, start=start(trimCoords), end=end(trimCoords))
    # Use IRanges coordinates to trim sequences and quality scores
    ids <- id(reads)
    qual <- SFastqQuality(qual) # reapply quality score type 
    bad <- which(width(seqs)==0)# looks for empty sequences
    
    ## removes empty sequences 
    if(length(bad)>0){
      seqs = seqs[-bad]
      qual = qual[-bad]
      ids = ids[-bad]
    }
    trimmed <- ShortReadQ(sread=seqs, quality=qual, id=ids)
    # Rebuilds reads object with trimmed sequences and quality scores
    # save, writes an external representation of R objects to the specified file
    save(file=paste(localFolder, "/", samplename, "-trimmed.RData", sep=""), list="seqs", "qual", "trimmed") 
    writeFastq(trimmed, file=paste( localFolder, "/", samplename, "-trimmed.fastq", sep=""))
    
    cat("Done.\n")
    logfile <- paste(logfile, "After Trimming: ", length(trimmed), " sequences \n", Sys.time(), "\n", sep ="")
    ## it can be useful to call gc after a large object has been removed, as this may prompt
    ## R to return memory to the operating system
    gc(verbose=F, reset=F)
    return.values = list(logfile, trimmed)
    return(return.values)
  }


########################
#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| Remove too short sequences      v1.02               |\n")
#cat("|                                                     |\n")
#cat("| (c) by Anja Lange, 2014                             |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("|                                                     |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")

remove_sequences_below_lengthcutoff <- function(data, sample.name, primer, with_primertable, min.length, logfile, localFolder){
  
  ## data are in ShortReadQ 
  minimum.length <- min.length
  
  if(with_primertable){### length of the primers has to be taken into account
    min.primer.length <- width(primer)
    minimum.length <- min.length + min.primer.length 
  }
  long.sequences <- data[width(data)>= minimum.length]
  
  writeFastq(long.sequences, file = paste(localFolder, "/", sample.name, "-wo_ShortSeq1.fastq", sep =""))
  
  removed.count <- length(data) - length(long.sequences)
  logfile <- paste(logfile, "After 1st length cutoff: ", length(long.sequences)," retained \n", Sys.time(), "\n", sep = "")
  cat( "removed  ", removed.count ," out of ", length(data), " sequences, that were below the length cutoff ", "\n", sep= "")
  gc(verbose=F, reset=F)
  return.value <- list(logfile, long.sequences)
  return (return.value) 
}
##########

#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| Deep Sequencing Quality checker v1.11               |\n")
#cat("|                                                     |\n")
#cat("| (c) by Dominik Heider, 2012                         |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("| dominik.heider@uni-due.de                           |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")

qualy <-
  function(data, cut_off_base, cut_off_seq, winsize, sample.name, logfile, localFolder){
    cat("Loading quality file...\n")
    
    reads <- data
    seqs <- sread(reads) # get sequence list
    lengthy <- length(seqs)
    qual <- PhredQuality(quality(quality(reads))) # get quality score list as PhredQuality
    ids <- id(reads)
    bad <- c() # remember bad quality sequences
    
    cat("Checking quality...\n")
    cat("[")
    list_counter <- c(1:100)
    .GlobalEnv[["list_counter"]] <- list_counter
    counter <- 1
    dots <- 100/length(seqs)
    .GlobalEnv[["counter"]] <- counter
    .GlobalEnv[["dots"]] <- dots
    
    ## parallel computing has to be executed from a shell and not a GUI !
    bad_list <- foreach(i = 1:length(qual)) %dopar% { 
      qualy_apply(i, qual, bad, cut_off_base, cut_off_seq, winsize) 
    }
    
    cat("]\n")
    cat("Saving files...\n")
    bad <- unlist(bad_list)
    if (length(bad)>0){
      seqs <- seqs[-bad,] # remove bad quality sequences
      qual <- qual[-bad,]
      ids <- ids[-bad,]
    }
    names(seqs) <- as.character(ids)
    
    outputfile = paste(localFolder, "/", sample.name, sep="")  
    ## writes a StringSet object to a FASTA file
    writeXStringSet(seqs, paste(outputfile, "-checked.fasta", sep = ""), append=TRUE)
    save(file=paste(outputfile, "-checked.RData", sep=""), list="seqs", "qual", "bad_list")
    
    cat(date())
    cat(" done.", "\n")
    
    logfile <- paste(logfile, "After Quality: ", length(seqs), " from initial ",
                     lengthy," sequences  kept. \n" , sep = "" )
    return.value<- list(logfile, seqs)
    gc(verbose=F, reset=F)
    return(return.value)
  }

#########################

qualy_apply <- 
  function(i, quality, bad_list, qualityCutoff, qualityCutoff.mean, winsize){
    count <- .GlobalEnv[["counter"]] ## integer starts with 1
    dots <- .GlobalEnv[["dots"]]  ##
    list_counter <- .GlobalEnv[["list_counter"]] ## list of intergers from 1 to 100
    if (!is.na(list_counter[1])){ ## if first value of list counter has a value
      if(floor(count) >= list_counter[1]){ ## rundet count ab
        cat(".")
        list_counter <- list_counter[-1] ## first element of list_counter is removed
      }
    }
    count <- count + dots ## adds value of dot, 2500 dots increase count of 1
    ## checks every quality score if it is lower than the cutt_off_base value returns boolean vector 
    ##with true where qaulity is low
    qualicheck <- as.integer(quality[i]) < qualityCutoff 
    ## calculates the mean of the quality scores
    meancheck <- mean(as.integer(quality[i]))
    
    if (winsize == FALSE){
      ## which, gives the TRUE indices of a logical object, allowing for array indices.
      ## contains the sequence bases with a quality below cut_off_base, or is the mean below
      ## cut_off_seq -> the index of the sequence is added to the bad_list
      if(length(which(qualicheck==TRUE)) > 0 | (meancheck < qualityCutoff.mean)){
        # adds index to the bad_list vector
        bad_list <- c(bad_list, i)
      }
    }
    else{
      qualvector <- as.vector(as.integer(quality[i]), mode="numeric")
      wincheck <- runmean(qualvector,winsize,endrule="mean",alg="C")
      if((wincheck < qualityCutoff.mean) | (length(which(qualicheck==TRUE)) > 0 )){
        bad_list <- c(bad_list, i)
      }
    }
    .GlobalEnv[["counter"]] <- count
    .GlobalEnv[["dots"]] <- dots
    .GlobalEnv[["list_counter"]] <- list_counter
    return (bad_list)
  }

###########################

primer_cutter <- 
  function(sequence, primer, max.length = NA){
    # sequence as DNAStringSet, primer as character
    
    match <- matchPattern(primer, sequence[[1]], fixed = FALSE) # returns XStringView
    if(length(match) == 0 || start(match) != 1){ ## did not match
      return(DNAStringSet()) ## return empty DNAStringSet
    }
    pureSeq.length <- width(sequence[1]) - width(primer)
    endPos <- end(match)
    if(pureSeq.length <= max.length){#the sequence is shorter than max.length
      cut.sequence <- subseq(sequence, start = endPos + 1, end = width(sequence)) ## DNAStringSet
    }else{
      cut.sequence <- subseq(sequence, endPos + 1, width = max.length) ## DNAStringSet
    }   
    return(cut.sequence)
  }

apply.DigDeeper(snakemake@input[["file"]], snakemake@input[["file_name"]], !!primer!!, !!conf_max_len!! ) <- function(file, sample.name, primer, max.length){
