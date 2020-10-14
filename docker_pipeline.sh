#!/bin/bash

if [ -z "$1" ]
then
	varname=$1
else
	varname=$PROJECT_NAME
fi
echo $varname
echo $1
echo $PROJECT_NAME
env_loc=$(conda info --base)/etc/profile.d/conda.sh

source $env_loc
conda activate snakemake
python demultiplexing.py "$varname"

varname+=".yaml"
cores=$(grep "cores : " $varname | awk '{print $3}')
python create_dataframe.py "$varname"

source $env_loc
conda activate snakemake
snakemake --use-conda --cores $cores --configfile $varname
exec sh
