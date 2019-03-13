#!/bin/bash

echo Please enter the name of the project

read varname
#source ~/anaconda3/etc/profile.d/conda.sh
conda activate snakemake
python demultiplexing.py "$varname"

varname+=".yaml"
cores=$(grep "cores : " $varname | awk '{print $3}')
python create_dataframe.py "$varname"

screen -S $varname bash -c "source ~/anaconda3/etc/profile.d/conda.sh;conda activate snakemake;snakemake --use-conda --cores $cores --configfile $varname; exec sh"
