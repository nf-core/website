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
Configuration of the pipeline is done via the `nextflow.config` file and any additional configuration files in `conf`. Typically, the `nextflow.config` contains all command-line parameters with their defaults, the handle for the docker container, nextflow options and the manifest, which informs about authors, nextflow and pipeline version and other things. The `conf` directory always contains a `base.config` file which is always loaded into `nextflow.config` and describes basic pipeline configurations, like CPU and memory usage for processes with low, medium and high requirements. Additionaly, `conf` contains a `igenomes.config` file which describes the locations of popular genomes that can be automatically downloaded for a pipeline run. Finally, the `test.config` and `test_full.config` files are test configurations that are loaded during test runs.

The `nextflow.config` defines different profiles that can be used to run the pipeline. These include the test profiles, which load, when used, the above mentioned configuration files into the `nextflow.config` file. Additionaly, the following profiles are available by default: `debug`, `conda`, `docker`, `singularity` and `podman`‚ which are used to run the pipeline in debug mode, using a conda environment or a using a container technology (docker, singularity or podman). 

## Continuous integration testing
To assure that nf-core pipelines don't break after some change is made to the code, we use automated continuous integration (CI) testing. This is done with GitHub actions, which are defined in the `.github/workflows` directory. Parameters and file paths are set in the `conf/test.config` and `conf/test_full.config`. Please see also HERE for how to set-up the test workflow for your pipeline.

## JSON Schema

## DSL2 and modules
Nextflow DSL2 allows for a more modularized design of pipelines and the reuse of components. Currently, most nf-core pipelines are still entirely written in DSL1, but in the near future all n-core pipelines will be written in DSL2. The nf-core team has developed a set of design patterns on how to best implement DSL2 pipelines, which should be used by all nf-core pipelines in order to assure standardization and the reuse of components. The following is meant to help understand certain design choices and how a nf-core DSL2 pipeline should be build.

### modules
Each pipeline has a `modules` directory which contains all the module code. A module here depicts a single process which involves - if possible - only a single tool/software. The `modules` directory is furthermore divided into `local`and `nf-core`sub-directories, which themselves each have a `process`and `subworkflow` directory. Modules contained in the `local` directory are specific to the pipeline, whereas `nf-core` modules are installed from the `nf-core/modules` repository. For instance, most pipelines that involve FastQ files will run the FastQC tool for quality control. The module needed for this can be easily reused from the `nf-core/modules` directory using the `nf-core/tools`package. The `process` directories contain modules which define single processes, which smaller workflows are contained in the `subworkflow` directories.

All modules furthermore load utility functions from a `functions.nf` script that must be contained in the `modules/local/process` directory. It contains simple functions to initialize the module options, get the software version, save files and get a path from a string. For further explanations of modules and how they should be structured in DSL2 pipelines see also HERE.

### The pipeline
The actual pipeline is now divided into two scripts, the `main.nf` script and the `<pipeline_name>.nf` script (e.g. `rnaseq.nf`). As with DSL1 pipelines, the script called when running the pipeline is the `main.nf` script. However, now all pipeline logic is stored in the `<pipeline_name>.nf` script, which includes all the necessary modules and subworkflows and executes them in the correct order. The bulk of code in this script will be to parse the command line parameters, load modules along with their specific parameters and execute them in the main workflow.

### Module parameters
One thing that might not be straightforward is how module parameters are handled in nf-core DSL2 pipelines. Every module and subworkflow, when loaded into the pipeline, has to be passed a groovy map containing module options. For single processes this is typically only a single `options` map, while for subworkflows these can be several maps that are then passed down to the correct processes within the subworkflows. These `option` maps are directly loaded from the `modules.config` file (contained in the pipeline `conf` directory), which is the place where all additional and optional parameters for modules are stored. Modules should be build in a way that they are flexible with respect to the parameters, so that most command line parameters can be passed to them via the `modules.config`. This way, all command line parameters and other options can be modified within a single script, which makes it easy for users to adjust the pipeline and at the same time makes modules resuable, as they stay flexible.

The `modules.config` file should contain a `params.modules` dictionary which lists every module used in the pipeline. For each module, the following fields can be specified:

- `args`: additional arguments appended to command in the module
- `args2`: Second set of arguments append to command in the module (multi-tool modules)
- `publish_dir`: Directory to publish the results
- `publish_by_id`: Publish results in separate folder by meta.id value
- `publish_files`: Groovy map where key = "file_ext" and value = "directory" to publish results for that file extension. The value of "directory" is appended to the standard "publish_dir" path as defined above. If publish_files == null (unspecified) all files are published. If publish_files == false no files are published.
- `suffix`: File name suffix for output files

### lib
Another feature of new nf-core pipelines, that is not strictly DSL2 related, is the `lib` directory. This directory simply contains several groovy-utility functions, that are used in the `main.nf` and and `<pipeline_name.nf>` scripts. 
