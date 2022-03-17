## All required dependencies are installed in the Gitpod environment.

## Print the command-line usage instructions for the *nf-core/viralrecon* pipeline

nextflow run nf-core/viralrecon --help


## In a new directory, run the *nf-core/viralrecon* pipeline with the provided test data

mkdir -p testrun
cd testrun
nextflow run nf-core/viralrecon -profile test,docker

## Try launching the Viralrecon pipeline using the *nf-core launch* command

nf-core launch viralrecon

## Download the *nf-core/viralrecon* pipeline for offline use using the *nf-core download* command

nf-core download viralrecon