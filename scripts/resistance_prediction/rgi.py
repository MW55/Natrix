import os
import subprocess
import sys

r1 = snakemake.input.r1
r2 = snakemake.input.r2
output_base = os.path.join("results", "rgi", f"{snakemake.wildcards.sample}_{snakemake.wildcards.unit}")
log_path = str(snakemake.log)

with open(log_path, "w") as logfile:
    base_cmd = [
        "rgi", "bwt",
        "-n", str(snakemake.params.threads),
        "-a", "bowtie2",
        "-o", output_base,
        "--local",
        "--clean"
    ]

    if r2 and len(r2) > 0:
        logfile.write("Paired-end detected. Running RGI bwt for paired-end data.\n")
        cmd = base_cmd + ["-1", r1, "-2", r2]
    else:
        logfile.write("Single-end detected. Running RGI bwt for single-end data.\n")
        cmd = base_cmd + ["-1", r1]

    logfile.write("Running command: " + " ".join(cmd) + "\n")
    try:
        subprocess.check_call(cmd, stdout=logfile, stderr=logfile)
    except subprocess.CalledProcessError as e:
        logfile.write("RGI command failed.\n")
        sys.exit(e)