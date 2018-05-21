#!/bin/bash

echo Please enter the name of the project

read varname
python demultiplexing.py "$varname"

varname+=".yaml"
cores=$(grep "cores : " $varname | awk '{print $3}')
python create_dataframe.py "$varname"

screen -S $varname bash -c "source activate snakemake;snakemake --use-conda --cores $cores --configfile $varname; exec sh"
