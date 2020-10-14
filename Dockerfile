FROM continuumio/miniconda3
COPY . /app
WORKDIR /app
RUN apt update && apt-get install -y libltdl7 && apt upgrade -y && apt-get purge -y && apt-get clean

# Create the environment:
#COPY snakemake.yaml .
RUN conda env create -f snakemake.yaml
RUN chmod +x docker_pipeline.sh
# Make RUN commands use the new environment:
#SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

#ENTRYPOINT ["bash"] 
#ENTRYPOINT ["./docker_pipeline.sh"]
CMD ["sh","-c", "./docker_pipeline.sh $PROJECT_NAME"]
