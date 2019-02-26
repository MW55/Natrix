from subprocess import call

if(len(snakemake.input)) == 2:
    output_edit = str(snakemake.output[0])[:-8]
    output_bad = str(snakemake.output[0])[:-8] + '_bad'
    call(['prinseq-lite.pl', '-verbose',
          '-fastq', snakemake.input[0], '-fastq2', snakemake.input[1], '-ns_max_n', '0',
          '-min_qual_mean', str(snakemake.params), '-out_good', output_edit,
          '-out_bad', output_bad, '-log', str(snakemake.log)])  #2>&1
else:
    output_edit = str(snakemake.output)[:-6]
    output_bad = str(snakemake.output)[:-6] + '_bad'
    call(['prinseq-lite.pl', '-verbose',
          '-fastq', snakemake.input[0], '-ns_max_n', '0', '-min_qual_mean', str(snakemake.params),
          '-out_good', output_edit, '-out_bad', output_bad, '-log', str(snakemake.log)])  #'2>&1'