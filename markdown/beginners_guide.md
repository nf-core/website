# What is nf-core?
A community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/docs/latest/index.html).
nf-core has three target audiences: facilities, single users and developers.
For facilities it provides highly automated and optimized pipelines that guaranty reproducibility of results for their users.
Single users profit from portable, documented and easy to use workflows.
But you can also become a developer and create your own pipeline using already available templates.

# What is Nextflow?
Nextflow is a *workflow manager*.
It has been developed specifically to ease the creation and execution of bioinformatics pipelines.
The benefits of having your pipeline in Nextflow include:

* Built-in GitHub support.
* Compatibility with most cluster executors.
* Docker support to eliminate the need for dependencies.
* Portability so you can run your pipeline anywhere laptop, cluster or cloud.

Whether your pipeline is a simple BLAST execution or a complex genome annotation pipeline, you can build it with Nextflow.

# What is Docker?
[Docker](https://www.docker.com/) is a tool to package and run applications locally in your machine within isolated environments called containers.
Importantly, a Docker container can be shared and, as such, all dependencies and versions required to run a specific script can be kept constant, even when a script is executed on a different machine, simply by running it in the respective Docker container.
The usage of [Docker in Nextflow](https://www.nextflow.io/docs/latest/docker.html) and container technology enables us to:

* Simplify the setup of complicated software, libraries and modules.
* Make our results, analysis and graphs 100% reproducible.
* Share our work with the world.

## When not to use Docker?
Docker is good when you are running your pipe locally. When you are working in a server it is recommended to work using [**Singularity**](https://www.sylabs.io/guides/2.5.1/user-guide/) containers.
Singularity is a container engine alternative to Docker. The main advantages of Singularity is that it can be used with unprivileged permissions, in this way you can run [Nextflow using Singularity](https://www.nextflow.io/docs/latest/singularity.html) in a server were you don't have root privileges.

## What is Conda?
[**Conda**](https://conda.io/) helps to find and install packages. More importantly, Conda easily creates, saves, loads and switches between environments. [Nextflow has built-in support for Conda](https://www.nextflow.io/docs/latest/conda.html) that allows run dependencies using Conda recipes and environment files. Also with this you can access and use the bioinformatic tools available in [Bioconda](https://bioconda.github.io/).

# How to run a pipeline?
In order to run a Nextflow pipeline from nf-core on your local computer you need to install Nextflow and Docker on your computer.

**System requirements:**
* Java 7 or 8
* Docker engine 1.10.x (or higher)
* recommended: 8GB of RAM


1. Install Nextflow

    `curl -s https://get.nextflow.io | bash`

2. Install nf-core tools

    `pip install nf-core`

3. You can check all the pipelines available by typing in your terminal

    `nf-core list`

4. To test that everything required is installed, run in your terminal

    `nextflow run nf-core/methylseq -profile test,docker -r dev`

5. Launch the pipeline of choice

    `nextflow run nf-core/<pipeline_name> --with-docker <parameters>`

For example, if you want to run `methylseq` pipeline, just type

`nextflow run nf-core/methylseq --reads 'path/*.fastq.gz' --outdir path/results --with-docker --genome < genome>`

You will find the specific parameters required for each pipeline in the documentation of the respective pipeline.
