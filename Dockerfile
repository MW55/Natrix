FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]

COPY . /app
WORKDIR /app
RUN apt update && apt-get install -y libltdl7 && apt upgrade -y && apt-get purge -y && apt-get clean

# Create the environment:
RUN conda env create -f natrix.yaml
RUN chmod +x docker_pipeline.sh
# Make RUN commands use the new environment:
#SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Creating envs for each task:
RUN env_loc=$(conda info --base)/etc/profile.d/conda.sh && source $env_loc && conda activate natrix && mkdir docker_dummy_env1 && touch docker_dummy_env1.csv && cp docker_dummyfiles/units.tsv docker_dummy.tsv && mkdir docker_dummy_env2 && touch docker_dummy_env2.csv && snakemake --configfile docker_dummyfiles/docker_dummy_env1.yaml --cores 1 --use-conda --conda-create-envs-only && snakemake --configfile docker_dummyfiles/docker_dummy_env2.yaml --cores 1 --use-conda --conda-create-envs-only && rm -rf docker_dummy_env1 && rm docker_dummy_env1.csv && rm -rf docker_dummy_env2 && rm docker_dummy_env2.csv && rm docker_dummy.tsv
CMD ["sh","-c", "./docker_pipeline.sh $PROJECT_NAME"]
