library('ShortRead')

# Check for primers & MIDs in the sequences
primer <- read.csv(snakemake@input[["primer_table"]], header=T, stringsAsFactors = F)
files_fwd <- list.files(snakemake@input[["data_folder"]], full.names=TRUE, pattern="\\R1.fastq.gz$")
files_rev <- list.files(snakemake@input[["data_folder"]], full.names=TRUE, pattern="\\R2.fastq.gz$")

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
text <- 'Starting MID/Primer checker\n'
cat(text, file=snakemake@log[[1]], append=T)
cat(text)
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
        cat(text, file=snakemake@log[[1]], append = TRUE)
        if (ignore == FALSE){
          cat(text)
          if(snakemake@config[['qc']][['all_primer']]){
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

#check if the reverse primer and/or the MID is still in the fastq sequences
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
          cat(text, file=snakemake@log[[1]], append = TRUE)
          if (ignore == FALSE){
            cat(text)
            if(snakemake@config[['qc']][['all_primer']]){
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
file.create('logs/qc_done')

