import subprocess
import sys

fq_files = snakemake.input.fq
index_file = snakemake.input.index
fq_prefix = snakemake.params.fq_prefix
bowtie2_params = snakemake.params.bowtie2_params
threads = snakemake.params.threads

# bowtie2 expects the base name of the index (without extension)
index_base = index_file
if index_file.endswith(".1.bt2"):
    index_base = index_file[:-6]

if len(fq_files) == 1:
    # Single-end: use the -U flag and --un for unaligned reads.
    output_file = fq_prefix + ".1"
    cmd = [
        "bowtie2",
        "--quiet",
        "-x", index_base,
        "-U", fq_files[0],
        "--un", output_file,
        "--threads", str(threads)
    ]
elif len(fq_files) == 2:
    # Paired-end: use the -1/-2 flags and --un-conc for unaligned pairs.
    cmd = [
        "bowtie2",
        "--quiet",
        "-x", index_base,
        "-1", fq_files[0],
        "-2", fq_files[1],
        "--un-conc", fq_prefix,
        "--threads", str(threads)
    ]
else:
    sys.exit("Error: Expected one (single-end) or two (paired-end) input fastq files, got {}.".format(len(fq_files)))

# If there are additional bowtie2 parameters, add them to the command.
if bowtie2_params:
    cmd.extend(bowtie2_params.split())

try:
    subprocess.check_call(cmd)
except subprocess.CalledProcessError as e:
    sys.exit("Error: bowtie2 command failed with exit code {}.".format(e.returncode))
