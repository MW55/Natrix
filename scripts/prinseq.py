from subprocess import call

# If the sequences are single-end, they still need the read identifier
# (R1), otherwise, the string slicing should remove the last six symbols
# instead of the last 8.
output_edit = str(snakemake.output[0])[:-8]
output_bad = str(snakemake.output[0])[:-8] + '_low_qual'

if(len(snakemake.input)) == 2:
    call(['prinseq-lite.pl', '-verbose',
          '-fastq', snakemake.input[0], '-fastq2', snakemake.input[1], '-ns_max_n', '0',
          '-min_qual_mean', str(snakemake.params), '-out_good', output_edit,
          '-out_bad', output_bad, '-log', str(snakemake.log)])
else:
    call(['prinseq-lite.pl', '-verbose',
          '-fastq', snakemake.input[0], '-ns_max_n', '0', '-min_qual_mean', str(snakemake.params),
          '-out_good', output_edit, '-out_bad', output_bad, '-log', str(snakemake.log)])