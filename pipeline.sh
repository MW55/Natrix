#!/bin/bash

echo Please enter the name of the project

read varname
env_loc=$(conda info --base)/etc/profile.d/conda.sh

source $env_loc
conda activate natrix

varname+=".yaml"
cores=$(grep "cores : " $varname | awk '{print $3}')
if python create_dataframe.py "$varname"; then
    screen -S $varname bash -c "source $env_loc;conda activate natrix;snakemake --use-conda --cores $cores --configfile $varname; exec sh"
fi
