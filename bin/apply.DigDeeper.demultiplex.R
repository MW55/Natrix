#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| apply.DigDeeper.demultiplex.R                       |\n")
#cat("|                                                     |\n")
#cat("| (c) by Anja Lange, 2014                             |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("| anja.lange@uni-due.de                               |\n")
#cat("| implements functions by Dominik Heider              |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")

apply.DigDeeper.demultiplex <- function(input_folder, type , max.length){
  
  input <- set_input_file(input_folder, type)
  
  filename <- paste(.GlobalEnv[["resultsFolder"]] , "/", input_folder, sep ="")
  
  if(with_primertable){
    ## creates a global variable with the primertable
    info <- read.csv(.GlobalEnv[["primerTable"]], header=T)
    .GlobalEnv[["info"]] <- info
  }
  
  if (type == "illumina"){
    trim_results <- trim(input, filename)
    a = date()
    cat("Trimming done.\n", file = .GlobalEnv[["logfile"]], append=T)
    cat(a, file = .GlobalEnv[["logfile"]], append=T)
    cat("\n\n", file = .GlobalEnv[["logfile"]], append=T)
  }else{
    trim_results <- input
  }
  
  cat("\n", "number of trimmed sequences: ", length(trim_results))
  
  if (min_length_of_sequence){
    first_removed <- remove_sequences_below_lengthcutoff(trim_results, info, with_primertable, min_length_of_sequence)
    a = date()
    cat("Short sequences removed (1).\n", file = .GlobalEnv[["logfile"]], append=T)
    cat(a, file = .GlobalEnv[["logfile"]], append=T)
    cat("\n\n", file = .GlobalEnv[["logfile"]], append=T)   
  }else{
    first_removed <- trim_results
  }
  
  .GlobalEnv[["first_removed"]]<- first_removed
  
  # Pakete der Größe 500.000 bauen, weil Illumina sonst den Speicher voll macht 
  # ceiling takes a single numeric argument x and returns a numeric vector containing
  # the smallest integers not less than the corresponding elements of x.
  num.packages = ceiling(length(first_removed) / 500000)
  for(i in 1:num.packages ){
    start = 500000 * i - 500000 + 1
    stop = 500000 * i
    if(stop>length(first_removed)){
      stop = length(first_removed)
    }
    cat(paste("Package number: ", i, "/", num.packages, "\n", sep=""))
    quality_results <- qualy(first_removed[start:stop], cut_off_base, cut_off_seq, winsize = FALSE , filename)
    print(start)
    print(stop)
  }
  a = date()
  cat("Quality checked.\n", file = .GlobalEnv[["logfile"]], append=T)
  cat(a,  file =.GlobalEnv[["logfile"]], append=T)
  cat("\n\n", file= .GlobalEnv[["logfile"]], append=T)
  
  quality_results <- readDNAStringSet(paste(filename, "-checked.fasta", sep=""))
  ##################
  
  if(length(quality_results)<2){
    stop("All sequences were bad in quality. Either change quality settings or consider new samples in lab.")
  }
  
  #quality_results are now in fasta format! Not in fastq
  if (with_primertable){
    ## removed contains all the sequences that were identified using the MID and primer, and these sequence parts are removed
    ## the sequence names correspond now to the sample
    removed <- removeForwardPrimer(quality_results, primerTable, filename, max.length)
    a = date()
    cat("Forward primer removed.\n", file=.GlobalEnv[["logfile"]], append=T)
    cat(a, file=.GlobalEnv[["logfile"]], append=T)
    cat("\n\n", file=.GlobalEnv[["logfile"]], append=T)  	
    cat("Forward primer removed.\n")
  
    ## test if sequences could be matched by the MID-primer sequences
    if(length(removed)<2){
      stop("All sequences were bad after primer remover. Consider new samples in lab.")
    }
    
    ## removes short sequences that are now below some threshold
    if (min_length_of_sequence){
      removed <- remove_short_ones_2(removed, filename, min_length_of_sequence)
      a = date()
      cat("Short sequences removed (2).\n", file=.GlobalEnv[["logfile"]], append=T)
      cat(a, file=.GlobalEnv[["logfile"]], append=T)
      cat("\n\n", file=.GlobalEnv[["logfile"]], append=T)
    }
    
    if(length(removed)<2){
      stop("All sequences were bad after length remover. Either change length settings or consider new samples in lab.")
    }
    ## creates fasta files for all samples and writes the sequences in the corresponding fasta files
    samples_names <- find_samples(removed, multicores, filename)
    
    a = date()
    cat("Samples deconvoluted.\n", file=.GlobalEnv[["logfile"]], append=T)
    cat(a, file=.GlobalEnv[["logfile"]], append=T)
    cat("\n\n", file=.GlobalEnv[["logfile"]], append=T)
  }	
  cat("All done!", "\n")
  }

#################################################################################################################################################
#############################      end dig_deeper    ############################################################################################
# data are stored in .GlobalEnv[["input"]] as ShortRead object

set_input_file <-  function(
  ### set the input file (known epitopes).
  path_to_file,
  ###f for 454 path to both files
  type_of_data
  ### the type of the data: 454 or Illumina
){
  ### known epitopes
  if (path_to_file != ""){
    if (type_of_data == "illumina"){
      ### illumina data are already in fastq format
      input <- read.illum(path_to_file) ###fastq-data
      filename <- strsplit(path_to_file, ".", fixed = TRUE)[[1]][1]
      .GlobalEnv[["filename"]] <- filename
    }
    else {
      input <- read.454(path_to_file) ###fastq-data
      filename <- strsplit(path_to_file, "/", fixed = TRUE)[[1]][2]
      .GlobalEnv[["filename"]] <- filename
    }
  }
  else {
    input <- c()
  }
  return(input)
}

### reads all FASTQ-formatted files in the directory path_to_file, returning a compact
### internal representation of the sequences and quality scores in the file
### return a ShortRead object
read.illum <- function(path_to_file){
  loaded_sequences <- readFastq(path_to_file, nrec=5)
  logfile <- .GlobalEnv[["logfile"]]
  cat("Starting with: ", length(loaded_sequences), "\n",file = logfile, append = TRUE)
  return (loaded_sequences)
}

### reads 454 files and returns fastq file
### 454 produces two seperate files, a fasta file with the sequences and a coresponding quality file
### informations from both files are combined in a fastq file
### needs a directory, not a file!!!
### if no path is provided a GUI for manual selection opens
read.454 <- function(path_to_directory = tcltk::tk_choose.dir()){
  rpath <- RochePath(experimentPath=path_to_directory)
  new_name <- gsub("/", "_", path_to_directory)
  reads <- read454(rpath) ###  read sequences and quality scores into a ShortReadQ instance.
  ### This class represents the directory location where Roche (454) result files (fasta sequences and qualities) can be found.
  ### read454 Pass arguments on to readFastaQual
  output <- paste(new_name, ".fastq", sep="")  #Concatenate vectors after converting to character.
  writeFastq(reads, file=output) 
  loaded_sequences <- readFastq(output, nrec=5)
  logfile <- .GlobalEnv[["logfile"]]
  cat("Starting with: ", length(loaded_sequences), "\n",file = logfile, append = T)
  return (loaded_sequences)
}

trim <- function(data, filename){
  #cat("\n")
  #cat("-------------------------------------------------------\n")
  #cat("| Poly-N trimmer v0.1                                 |\n")
  #cat("|                                                     |\n")
  #cat("| (c) by Dominik Heider, 2012                         |\n")
  #cat("| University of Duisburg-Essen, Germany               |\n")
  #cat("| dominik.heider@uni-due.de                           |\n")
  #cat("|                                                     |\n")
  #cat("-------------------------------------------------------\n")
  #cat("\n")
  cat("Loading Fastq file...\n")
  reads <- data
  logfile <- .GlobalEnv[["logfile"]]
  cat("Before Trimming: ", length(reads), "\n",file = logfile, append=TRUE)
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
  save(file=paste(filename, "-trimmed.RData", sep=""), list="seqs", "qual", "trimmed") 
  writeFastq(trimmed, file=paste(filename, "-trimmed.fastq", sep=""))
  cat("Done.\n")
  cat("After Trimming: ", length(trimmed), "\n",file = logfile, append = TRUE)
  ## it can be useful to call gc after a large object has been removed, as this may prompt
  ## R to return memory to the operating system
  gc(verbose=F, reset=F)
  return (trimmed)
}

remove_sequences_below_lengthcutoff<- function(data, primer.table, with_primertable, min.length ){
  
  #cat("\n")
  #cat("-------------------------------------------------------\n")
  #cat("| Remove too short sequences      v1.01               |\n")
  #cat("|                                                     |\n")
  #cat("| (c) by Anja Lange, 2014                             |\n")
  #cat("| University of Duisburg-Essen, Germany               |\n")
  #cat("|                                                     |\n")
  #cat("|                                                     |\n")
  #cat("-------------------------------------------------------\n")
  #cat("\n")
  ## data are in ShortReadQ 
  minimum.length <- min.length
  
  if(with_primertable){### length of the primers has to be taken into account
    ### combines poly_N, MID and primer to one String, for all samples
    primer.sequences <- sapply(X = seq(1:length(primer.table[,1])), FUN = function(x){
      paste(  primer.table[x, "poly_N"] , primer.table[x, "MID"], primer.table[x, "specific_forward_primer"], sep ="")
    })
    ## calculates the length of all primer sequences including MID and poly_N
    primer.length <- sapply(X = seq(1:length(primer.sequences)), FUN = function (y){
      prim <- unlist(strsplit(primer.sequences[y], split = ""))
      prim.length <- length(prim)
      return(prim.length)
    })
    min.primer.length <- min(primer.length)
    minimum.length <- min.length + min.primer.length 
  }
  long.sequences <- data[width(data)>= minimum.length]
  
  writeFastq(long.sequences, file=paste(.GlobalEnv[["resultsFolder"]],"/-removed_short_Sequences.fastq", sep =""))
  
  removed.count <- length(data) - length(long.sequences)
  cat("\n", "After 1st length cutoff: ", length(long.sequences),"\n", file = .GlobalEnv[["logfile"]], append = TRUE)
  cat("\n", "removed  ", removed.count ," out of ", length(data), " sequences, that were below the length cutoff", "\n")
  gc(verbose=F, reset=F)
  return (long.sequences) 
}


qualy <- function(data, cut_off_base, cut_off_seq, winsize, filename){
  #cat("\n")
  #cat("-------------------------------------------------------\n")
  #cat("| Deep Sequencing Quality checker v1.11               |\n")
  #cat("|                                                     |\n")
  #cat("| (c) by Dominik Heider, 2012                         |\n")
  #cat("| University of Duisburg-Essen, Germany               |\n")
  #cat("| dominik.heider@uni-due.de                           |\n")
  #cat("|                                                     |\n")
  #cat("-------------------------------------------------------\n")
  cat("\n")
  cat("Loading quality file...\n")
  
  outputfile = paste(filename, "-checked.fasta", sep="")  
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
  bad_list <- foreach(i = 1:length(qual)) %dopar% { qualy_apply(i, qual, bad, cut_off_base, cut_off_seq, winsize) }
  cat("]\n")
  cat("Saving files...\n")
  bad <- unlist(bad_list)
  if (length(bad)>0){
    seqs <- seqs[-bad,] # remove bad quality sequences
    qual <- qual[-bad,]
    ids <- ids[-bad,]
  }
  names(seqs) <- as.character(ids)
  
  ## writes a StringSet object to a FASTA file
  writeXStringSet(seqs, outputfile, append=TRUE)
  save(file=paste(outputfile, ".RData", sep=""), list="seqs", "qual", "bad_list")
  done <- date()
  cat(done)
  cat("Done.\n")
  logfile <- .GlobalEnv[["logfile"]]
  cat("After Quality: ", length(seqs), "\n",file = logfile, append = TRUE)
  cat("From initial ", lengthy, length(seqs), " kept.\n")
  gc(verbose=F, reset=F)
  return(seqs)
}

qualy_apply <-function(i, quality, bad_list, qualityCutoff, qualityCutoff.mean, winsize){
  count <- .GlobalEnv[["counter"]] ## integer starts with 1
  dots <- .GlobalEnv[["dots"]]  ##2e-04
  list_counter <- .GlobalEnv[["list_counter"]] ## list of intergers from 1 to 100
  if (!is.na(list_counter[1])){ ## if first value of list counter has a value
    if(floor(count) >= list_counter[1]){ ## rundet count ab
      cat(".")
      list_counter <- list_counter[-1] ## first element of list_counter is removed
    }
  }
  count <- count + dots ## adds value of dot, 5000 dots increase count of 1
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

removeForwardPrimer<-function(data, primerTable, filename,  max.length){
  #cat("\n")
  #cat("-------------------------------------------------------\n")
  #cat("| Primer Remover (F) v1.1                             |\n")
  #cat("|                                                     |\n")
  #cat("| (c) by Dominik Heider, 2012                         |\n")
  #cat("| University of Duisburg-Essen, Germany               |\n")
  #cat("| dominik.heider@uni-due.de                           |\n")
  #cat("| optimized by Anja Lange 2014                        |\n")
  #cat("-------------------------------------------------------\n")
  #cat("\n")
  cat("Loading Fasta file...\n")
  info <- .GlobalEnv[["info"]]
  sequences <- data
  .GlobalEnv[["sequences"]] <- data
  comb = c()
  ## concatenates the MID and forward-primer sequences for sample/line in the table  
  for(i in 1:length(info[,1])){
    comb <- c(comb, paste(info[i, "poly_N"], info[i,"MID"], info[i,"specific_forward_primer"], sep=""))
  }
  info <- cbind(info, comb)
  cat("Removing primer...\n")
  .GlobalEnv[["data.new"]] = c() 
  
  cat("[")
  ### Packages each with 10,000 sequences:
  number_of_sequences <- length(data)
  packages <- c()
  if (number_of_sequences <= 10000){
    packages <- c(1,number_of_sequences)
  }
  else {
    pieces <- ceiling(number_of_sequences/10000)
    first_step <- 1
    for (i in 1:pieces){
      next_step <- first_step + 9999
      if (next_step > number_of_sequences){
        next_step <- number_of_sequences
      }
      ## integer vector, every two integers specify the range of one package
      packages <- c(packages, first_step, next_step)
      first_step <- next_step + 1
    }
  }
  
  ## iterates over all packages of size 10000
  ## all_data will contain a list with the shortend sequences and a DNAStringSet with the sequences
  ## that were not found
  ## .combine specifies the function that is used to process the tasks results as they generated; here
  ## from all iterations the new_data DNAStringSets are appended, contain the identified and shortened sequences
  ## and the sequence-packages that contain those sequences that couldn't be matched by a primer
  cat("\n")
  cat("start iterate over packages of 10000: ", date())
  cat("\n")
  
  package.number <- length(packages)/2
  package.counter <- 1
  all_data <- foreach (p = seq(1, length(packages)-1, by=2), .combine="last_one") %do% {
    cat("\n", "package ", package.counter ,"/", package.number )
    cat("\n", date())
    package.counter <- package.counter + 1
    pa <- data[packages[p]:packages[p+1]]   ## sequences of the package
    .GlobalEnv[["sequences_package"]] <- pa ## package
    .GlobalEnv[["data.new"]] <- c()         ## empty vector
    ## iterates over all samples, lines in primer table
    ## the sequence packages get smaller during the iteration, matches sequences are removed
    foreach(i = 1:length(comb)) %do% { remove_forward_primer_apply(i, comb, info, max.length) }
    sequences_package <- .GlobalEnv[["sequences_package"]] #
    new_data <- .GlobalEnv[["data.new"]]
    ## returns list with new_data and sequence_package and applies the function "last_one"to all returned list
    ## in the for each loop
    return (list(new_data, sequences_package)) 
  }
  cat("]\n")
  cat("end iterate over packages of 10000: ", date())
  cat("\n")
  ## 
  sequences_worked <- all_data[[2]] ## weird name, contains sequences that didn't match a primer
  new_data <- all_data[[1]] ## identified and shortened sequences
  writeXStringSet(sequences_worked, paste(filename, "-bad_primer.fasta", sep="")) #, open = "w", nbchar = 60)
  cat(paste(length(sequences_worked), " sequences of ", (length(new_data)+length(sequences_worked)), " sequences removed.\n", sep=""))
  ## stores matched and shortened data in the variable test
  test <- new_data
  rm(new_data)
  writeXStringSet(test, paste(filename,"-prim_removed.fasta", sep=""))#, open = "w", nbchar = 60)
  cat("Done.\n")
  logfile <- .GlobalEnv[["logfile"]]
  cat("After Primer Remover: ", length(test), "\n",file = logfile, append = TRUE)  
  gc(verbose=F, reset=F)
  return(test)
}

last_one <- function(accum, ...){
  list(append(accum[[1]], ...[[1]]), append(accum[[2]], ...[[2]]))
}

remove_forward_primer_apply <- function(i, comb, info, max.length){
  ## comb contains character vector with poly_N , MID and forward primer sequence for each sample
  ## info contains the primer table, sequence_package, and seq are DNAStringSet Objects of size 10000
  ## i specifies the the sample/ MID-Primer combination the function is looking for 
  
  seq <- .GlobalEnv[["sequences_package"]] #10000 sequences, for the first primer 
  new_data <- .GlobalEnv[["data.new"]] ## for i=1 an empty vector, then appends seq[hits] with removed primers
  
  if (length(seq)>0){
    ##------------ changed this part---------
    ## looks for poly-N + MID+ Primer, poly N is,an IUPAC ambiguity code in the pattern
    ## can match any letter in the subject that is associated with the code
    ##!!problem if sequences is of low quality and has a lot of Ns
    hits <- vmatchPattern(comb[i], seq, fixed = FALSE) ## returns a ByPos_MIndex object
    hits.count <- vcountPattern(comb[i], seq, fixed = FALSE) ##numeric vector with 1 for a match, else 0
    ### calculates a vector with the indices of the sequences that were matched, 0 for not matched sequences
    hits.index <- sapply(X= seq(1:length(hits.count)), FUN= function(k){
      if(hits.count[k]==1){
        value <- k
      }else{
        value <-0
      }
      return(value)
    })
    hits.index <-hits.index[hits.index != 0] ### numeric vector contains only indices of positive matches
    ##------------------------------------------  
    if(!is.na(hits.index[[1]])){
      .GlobalEnv[["seq"]] <- seq ## all sequences of the package/ DNAStringSet
      
      sample.name <- as.character(info[i, "Probe_FWD"])
      #iterates over all hits, cuts off the matched sequences 
      ## shortens the sequence in the .GlobalEnv[["seq"]] and changes the name of the sequence 
      foreach(j = 1:length(hits.index)) %do% { remove_forward_primer_cut_matched_pattern(i, hits.index[j], hits[[hits.index[j]]], sample.name, max.length) }
      ## global variable seq contains shortened sequence
      seq <- .GlobalEnv[["seq"]]
      ## all sequences where the primers are removed are appended to new_data
      new_data <- append(new_data, seq[hits.index])
      #sequences with no hits are kept in sequence_package, is used to look for other primers
      .GlobalEnv[["sequences_package"]] <- seq[- hits.index]
    }
  }
  .GlobalEnv[["data.new"]] <- new_data
  
  gc(verbose=F, reset=F)
  return (new_data)
}

#cat("\n")
#cat("-------------------------------------------------------\n")
#cat("| remove_forward_primer_cut_matched_pattern (F) v0.1  |\n")
#cat("|                                                     |\n")
#cat("| (c) by Anja Lange    , Jan 2014                     |\n")
#cat("| University of Duisburg-Essen, Germany               |\n")
#cat("|                                                     |\n")
#cat("|                                                     |\n")
#cat("-------------------------------------------------------\n")
#cat("\n")
remove_forward_primer_cut_matched_pattern<- function(primerIndex, hitIndex, matchedPos, sample.name, max.length){
  ## primerIndex, numeric value, position of the matched primer in the primer list
  ## hitIndex, numeric value , gives the position of the matched sequences in the sequence package
  ## matchedPos: IRanges object, contains the precise matched position in the sequence
  ## sequences: DNAStringSet object with all sequences
  ## max.length: specified maximal sequence length
  index <- hitIndex
  matched.sequence <- .GlobalEnv[["seq"]][index] ##DNAStringSet
  seq.name <- names(matched.sequence) 
  new.name <- sample.name
  
  cutpoint <- (end(matchedPos)+1)
  seq.length <- width(matched.sequence)
  final.length <- seq.length - end(matchedPos)
  ## checks if the sequence would be longer then the specified max.lenght
  if(final.length <= max.length){
    cut.sequence <- subseq(matched.sequence, start = cutpoint, end = seq.length ) ## DNAStringSet
  }else{
    cut.sequence <- subseq(matched.sequence, start= cutpoint, width = max.length ) ## DNAStringSet
  }
 
  ## replaces inserts the cut sequence in the sequence list
  .GlobalEnv[["seq"]][index] <- cut.sequence
  names(.GlobalEnv[["seq"]])[index]<- new.name
}


remove_short_ones_2 <- function(quality_results, filename, min_length_of_sequence){
  #cat("\n")
  #cat("-------------------------------------------------------\n")
  #cat("| Remover (short Seqs) 2 v1.11                        |\n")
  #cat("|                                                     |\n")
  #cat("| (c) by Bettina Budeus, 2012                         |\n")
  #cat("| University of Duisburg-Essen, Germany               |\n")
  #cat("| bettina.budeus@stud.uni-due.de                      |\n")
  #cat("|                                                     |\n")
  #cat("-------------------------------------------------------\n")
  #cat("\n")
  #cat("Loading quality file...\n")
  #new_filename <- gsub("[ -]", "_",strsplit(filename, ".fasta"))
  #outputfile = paste(new_filename, "_shortend_2.fasta", sep="")
  outputfile = paste(filename, "_shortend_2.fasta", sep = "")
  result <- (quality_results[width(quality_results)>=min_length_of_sequence])
  writeXStringSet(result, outputfile)
  cat("Done.\n")
  logfile <- .GlobalEnv[["logfile"]]
  cat("After shorter sequence remover 2: ", length(result), "\n",file = logfile, append = TRUE)
  cat("From : ", length(quality_results), length(result), "kept.")
  return (result)
}


find_samples <- function(data, multicores, filenames){
  #cat("\n")
  #cat("-------------------------------------------------------\n")
  #cat("| Find Samples (F) v2.0                               |\n")
  #cat("|                                                     |\n")
  #cat("| (c) by Dominik Heider, 2013                         |\n")
  #cat("| University of Duisburg-Essen, Germany               |\n")
  #cat("| dominik.heider@uni-due.de                           |\n")
  #cat("|                                                     |\n")
  #cat("-------------------------------------------------------\n")
  #cat("\n")
  logfile <- .GlobalEnv[["logfile"]]
  cat("Now generating single files", "\n",file = logfile, append = TRUE)
  cat("Loading sequence...\n")
  sequences <- data
  namen <- names(sequences) ##string vector,  extracts for all sequences the names/identifiers
  dateien = unique(namen) ### removes all duplicate names
  cat(paste(" identified ", length(dateien), " samples: \n", sep = ""),file = logfile, append = TRUE)
  cat(paste0(dateien, collapse = ", "), "\n", file = logfile, append = TRUE)
  
  ## creates for each sample a seperate resulte/folder
  foreach(i = 1:length(dateien)) %do% {
    tmp = dateien[i]
    local.resultsFolder = paste(.GlobalEnv[["resultsFolder"]], "/", tmp, sep = "")
    sys = system(paste("mkdir -p ./", local.resultsFolder, sep = "" ))
  }  
  ## iterates over all sample names
  all_names <- foreach(i = 1:length(dateien)) %dopar% { find_samples_apply(sequences, dateien[i], filenames) }
  gc(verbose=F, reset=F)
}

## i indexes the name that is looked for
## sequences: DNAStringSet with all good sequences
## namen: sequences names without duplicates
## filenames 

find_samples_apply <- function(sequences, sample.name, filenames){
  local.resultsFolder = paste(.GlobalEnv[["resultsFolder"]], "/", sample.name, sep = "")
  ## filename contains the directory path and name of the new file
  filename = paste(local.resultsFolder, "/", sample.name, ".fasta", sep="")
  
  ## vector with names of all sequences
  tmp = names(sequences)
  ## selects those that have the same name as the indexed name
  ## returns a vector with the true indices
  hits = which(tmp == sample.name)
  
  ## selects only the sequences with indexed name
  sequenzen <- sequences[hits]
  name <- names(sequenzen)
  
  ## a file will be reported as existing only if you have the permissions needed by stat
  ### writes the found sequences, with the corresponding sample name in a fasta file 
  if(file.exists(filename)){
    writeXStringSet(sequenzen, filename, append=TRUE)
  }else{
    writeXStringSet(sequenzen, filename)
  }
  return(name)
}
