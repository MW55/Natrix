system("mkdir -p results/assembly")
results_dir <- "results/assembly"
#step 1 #####
read.Assembler <- function(input_folder, threshold, minoverlap, minqual, minlength, maxlength, primer.offset, pre.filter, pre.filter.mq){

  ## uses pandaseq to assemble paired End reads  
  sys = system("pandaseq -v 2>&1 ", intern = T)
  cat(sys[1], "\n\n", file = snakemake@log[[1]] , append = T)

   
    cat(paste(Sys.time(), "| Start assembling ", file.name, "\n", sep = ""), file = snakemake@log[[1]] , append = T)
    cat("read.Assembler: \n")
    if(file.name == file.name.r){
      this.primer <- primer[primer$Probe_FWD == file.name,]
      this.primer.f <- as.character(this.primer$f_primer)
      this.primer.r <- as.character(this.primer$r_primer)
      
      panda.out <- paste(localFolder, "/", file.name, "_assembled.fastq", sep ="")
      
      cat("Checking for existence of some files in ", localFolder, ": \n")
      cat("Forward: ", file.exists(file.f), " ",file.f, "\n")
      cat("Reverse: ", file.exists(file.r), " ", file.r, "\n") 
      
      ## call PandaSeq
      if(snakemake@config[["qc"]][["all_primer"]]){
        cmd = paste("pandaseq -f ", file.f, " -r ", file.r, " -B -a -F -g ", snakemake@log[[2]], " -w ", panda.out, " -N  -t ", threshold, " -o ", minoverlap, " -l ", minlength,
                                        " -L ", maxlength, " -C min_phred:", minqual , sep= "" )
      }else{
        cmd = paste("pandaseq -f ", file.f, " -r ", file.r, " -B -a -F -g ", snakemake@log[[2]], " -w ", panda.out, " -N -p \"",
                                        this.primer.f, "\" -q \"", this.primer.r, "\" -t ", threshold, " -o ", minoverlap, " -l ", minlength,
                                                            " -L ", maxlength, " -C min_phred:", minqual , sep= "" )
      }
      
      sys = system(cmd, intern = T) 
      cat(paste(Sys.time(), "| Statement executed: \n", cmd, "\n\n", sep=""), file = snakemake@log[[1]], append = T)
      
      cat("panda logging: ", file.exists(snakemake@log[[2]])," ", snakemake@log[[2]], "\n")
      cat("panda assembled: ", file.exists(panda.out), " ", panda.out, "\n\n")
      
      ### read fastq file, rename sequences and save as fasta
      if(file.info(panda.out)$size > 0){
        assembled.fastq <- readFastq(panda.out)
        assembled.seq <- sread(assembled.fastq)
        
        nr <- length(assembled.fastq)
        cat(paste(Sys.time(), "| assembled ", nr, " reads \n \n", sep=""), file = snakemake@log[[1]] , append = T)
        new.names <- sapply(X = 1: nr, FUN = function(x, file.name){
          name <- paste(file.name, "_", x, sep= "")
          return(name)
        }, file.name = file.name)
        
        names(assembled.seq) <- new.names
        writeXStringSet(assembled.seq, filepath = paste(localFolder, "/", file.name , ".fasta", sep=""))

      }else{ ## no reads assembled
       # main.folder <- strsplit(snakemake@output[["dir"]], split = "/results")[1] 
        ## make folder for sample with no assembled reads
        main.folder <- '/results' 
        bad.assembled.folder <- paste(main.folder, "/not.assembled", sep ="")
        if(! file.exists(bad.assembled.folder)){
          cmd = paste("mkdir -p ", bad.assembled.folder, sep ="" )
          sys = system(cmd, intern = T)   
        }        
        ## move folder
        cmd = paste("mv", localFolder,  bad.assembled.folder, sep =" ")
        sys = system(cmd, intern = T)
        message = "no reads assembled \n\n" 
        cat(message, file = snakemake@log[[1]] , append = T)
        
      }
      
    }else{
      stop("forward and reverse file don't match!")
    }     
  } 
}





if(snakemake@config[["merge"]][["paired_End"]]){
#### assembly of paired end reads
#### function can be found in fun.R
read.Assembler(snakemake@input[["data_folder"]], snakemake@config[["qc"]][["threshold"]], snakemake@config[["qc"]][["minoverlap"]], snakemake@config[["qc"]][["minqual"]], snakemake@config[["qc"]][["minlen"]], snakemake@config[["qc"]][["maxlen"]], snakemake@config[["qc"]][["primer_offset"]], snakemake@config[["qc"]][["prefilter"]], snakemake@config[["qc"]][["mq"]])
}else{ ## singleEnd
source("./bin/apply.DigDeeper1.03.R")
## load primer table
primer  <- read.csv(snakemake@input[["primer_table"]], header=T, stringsAsFactors = F)
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
  
  apply.DigDeeper(file, file.name, this.primer, snakemake@config[["qc"]][["maxlen"]])
}
rm(files)
}
gc()
file.create('logs/assembly_done')
