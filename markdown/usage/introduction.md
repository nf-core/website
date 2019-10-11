# WE NEED A TOP LEVEL HEADER HERE, AND IT NEEDS TO BE THE FIRST LINE

## What is nf-core

nf-core is a community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/docs/latest/index.html).

nf-core has three target audiences: facilities, single users and developers.
For facilities it provides highly automated and optimized pipelines that guaranty reproducibility of results for their users.
Single users profit from portable, documented and easy to use workflows.
But you can also become a developer and write your own pipeline in Nextflow using already available templates and helper tools.

## What is Nextflow

Nextflow is a *workflow manager*.
It has been developed specifically to ease the creation and execution of bioinformatics pipelines.
The benefits of having your pipeline in Nextflow include:

* Built-in GitHub support.
* Compatibility with virtually all computational infrastructures, including all major cluster job schedulers.
* Integrated software dependency management (Docker, Singularity, Conda).
* Portability so you can run your pipeline anywhere: laptop, cluster or cloud.
* Reproducibility of analyses independent of time and computing platform.

Whether your pipeline is a simple BLAST execution or a complex genome annotation pipeline, you can build it with Nextflow.

To learn more about nextflow, we have adapted a nextflow tutorial which you can try: [nextflow tutorial](/usage/nextflow_tutorial).

## Software requirements

An analysis pipeline chains the execution of multiple tools together.
In order for this to work, each of those software tools must be installed.
Historically, this can be a source of great frustration and a key step where reproducibility between analyses is lost.
All nf-core pipelines fully utilise the built-in support for software packaging that Nextflow offers.

### Docker

[Docker](https://www.docker.com/) is a tool to package and run software within isolated environments called containers.
Importantly, a Docker container images can be shared. As such, all dependencies and versions required to run a specific script can be kept constant, even when a script is executed on a different machine, simply by running it in the respective Docker container.
The usage of [Docker in Nextflow](https://www.nextflow.io/docs/latest/docker.html) and container technology enables us to:

* Simplify the setup of complicated software, libraries and modules.
* Make our results, analysis and graphs 100% reproducible.
* Share our work with the world.

### Singularity

It's not always possible to run Docker - it requires special system permissions which may not be feasible in a shared computing environment.
[Singularity](https://www.sylabs.io/guides/3.1/user-guide/) is a container engine alternative to Docker designed to run large analysis jobs on high performance compute clusters.
The main advantage is that it can be used with unprivileged permissions: in this way you can run [Nextflow using Singularity](https://www.nextflow.io/docs/latest/singularity.html) on a server were you don't have root privileges.

As with docker, singularity allows us to bundle all software requirements together into a single image which comes with each nf-core pipeline. This gives simplicity and reproducibility.

### Conda

[Conda](https://conda.io/) is a software packaging system that helps to find and install packages.
Conda can create isolated 'environments' with software installations that can be switched between.

[Nextflow has built-in support for Conda](https://www.nextflow.io/docs/latest/conda.html) that allows us to manage dependencies using Conda recipes and environment files. The nf-core pipelines make extensive use of the bioinformatics conda channel [Bioconda](https://bioconda.github.io/) making it simple to install all required tools.

Conda makes the installation and management of software dependencies far simpler than manual installations. It ensures correct version pinning and so represents good reproducibility.
However, the software still runs in your operating system environment. As such, if it's possible to use either docker or singularity, those options are preferred.

## How to run a pipeline

In order to run a Nextflow pipeline from nf-core on your local computer you need to install Nextflow and Docker on your computer.

System requirements:

* Java 7 or 8
* _recommended:_ at least 8GB of RAM
* _recommended:_ Docker engine 1.10.x (or higher)

1. Install Nextflow

    ```bash
    curl -s https://get.nextflow.io | bash
    ```

2. Install nf-core tools

    ```bash
    pip install nf-core
    ```

3. You can check all the pipelines available by typing in your terminal

    ```bash
    nf-core list
    ```

4. To test that everything required is installed, try running a pipeline test in your terminal

    ```bash
    nextflow run nf-core/methylseq -profile test,docker
    ```

5. Launch the pipeline of choice. `parameters` are optional, when they are not specified they are taken from the `nextflow.config` file.

    ```bash
    nextflow run nf-core/<pipeline_name> -profile standard,docker [parameters]
    ```

    For example, if you want to run `methylseq` pipeline, just type:

    ```bash
    nextflow run nf-core/methylseq -profile standard,docker -reads 'path/*.fastq.gz' --outdir path/results --genome <genome>
    ```

    You will find the specific parameters required for each pipeline in the documentation of the respective pipeline.

## Memory usage

We recommend to check your possible memory usage limit on a machine where you are going to launch a pipeline, as a Java machine used to orchestrate the pipeline will require a certain amount of physical memory. For example, if resources on the machine are managed with [cgroups](http://man7.org/linux/man-pages/man7/cgroups.7.html), you could check the amount of available memory using the following commands:

```bash
systemctl -t slice | grep type-your-username-here
```

This will give you your username in `cgroups`, for example:

```bash
user-30030283.slice                           loaded active active User Slice of username
```

Next, use your `cgroups` username (in this example, `user-30030283.slice`) to inquire about a possible memory usage limit:

```bash
systemctl status user-30030283.slice | grep "Memory:"
```

It will give you something like this:

```bash
   Memory: 24.2M (limit: 8.0G)
             ├─ 2816 grep --color=auto Memory:
```

So, the amount of memory available to this user on this machine is 8 GB.

> NB. We recommend to have at least 4 GB available.

You may also try to explicitly control the Java heap size by setting `NXF_OPTS='-Xms1g -Xmx4g'` (`1g` and `4g` are just an example!) in your command line environment (typically, `~/.bashrc` or `~./bash_profile`).

> NB. Your values for the `-Xms` and `-Xmx` arguments should be less then the amount of memory available to you.

See [-Xms](https://docs.oracle.com/cd/E13150_01/jrockit_jvm/jrockit/jrdocs/refman/optionX.html#wp999527) and [-Xmx](https://docs.oracle.com/cd/E13150_01/jrockit_jvm/jrockit/jrdocs/refman/optionX.html#wp999528) Java machine arguments for details.

## Helper tools

To help you manage your nf-core pipelines and discover updates, we have written some command-line helper tools.
These allow you to list all available pipelines and versions, with information about what versions you're running locally.
There are also commands to help downloading pipelines for use offline.

To find out more about these tools, read the [Tools](/tools) page.
