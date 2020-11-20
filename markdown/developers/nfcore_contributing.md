# nf-core contribution guidelines

Nf-core pipelines adhere to common design principles that have been crafted through discussions and extensive real-world experiences with pipelines over time. We understand, however, that not all about nf-core pipelines is obvious and can be confusing to people new to nf-core. To help you along getting started with nf-core, writing your pipeline or contribute otherwise, we have created this documentation. Hopefully you will find most things that might confuse you described here, and if not, be sure to bug the nf-core Slack channel!

## General stuff
The best way to write a new pipeline is to start from the template (see here). After doing that, you will have a number of files in your repository with different purposes. Here's what they are for:

### Files
* `main.nf`: This is the main nextflow file which will get executed if the pipeline is run. It is the most important file and contains the pipeline description.

* `nextflow.config`: The main nextflow configuration file. It contains the default pipeline parameters, nextflow configuration options and information like pipeline and nextflow version, among others.

* `README.md`: Basic information about the pipeline and usage

* `environment.yml`: Holds the conda environment specification. This is used for building the conda environment if the pipeline is used with conda but also for building the docker image.

* `Dockerfile`: Instrunctions for building the pipeline docker image. If possible, this will just import the `environment.yml` file and build a conda environment in the container.

* `nextflow_json.schema`: Hell if I know what that's for ...

* `CHANGELOG.md`: The change log contains information about, well, changes made to the pipeline.

* `LICENSE`: The license - should be MIT

* `CODE_OF_CONDUCT.md`: The nf-core code of conduct.

* .gitignore: Lists GitHub settings

* .gitattributes: Lists files that should be ignored by GitHub.

### Directories

* `.github`: Other GitHub specific files, e.g. for specifying templates and GitHub actions
* assets: Any additional files needed for the pipeline
* bin: Contains scripts to be executed in the pipeline
* conf: Additional configuration profiles that can be loaded in the `nextflow.config`file
* docs: Markdown files for documenting the pipeline

## Pipeline configuration
Configuration of the pipeline is done via the *nextflow.config* file and any additional configuration files in *conf*. Typically, the *nextflow.config* contains all command-line parameters with their defaults, the handle for the docker container, nextflow options and the manifest, which informs about authors, nextflow and pipeline version and other things. The *conf* directory always contains a *base.config* file which is always loaded into *nextflow.config* and describes basic pipeline configurations, like CPU and memory usage for processes with low, medium and high requirements. Additionaly, *conf* contains a *igenomes.config* file which describes the locations of popular genomes that can be automatically downloaded for a pipeline run. Finally, the *test.config* and *test_full.config* files are test configurations that are loaded during test runs.

The *nextflow.config* defines different profiles that can be used to run the pipeline. These include the test profiles, which load, when used, the above mentioned configuration files into the *nextflow.config* file. Additionaly, the following profiles are available by default: *debug*, *conda*, *docker*, *singularity* and *podman* which are used to run the pipeline in debug mode, using a conda environment or a using a container technology (docker, singularity or podman). 

## Continuous integration testing

## Configuration



## JSON Schema


## DSL2 and modules

