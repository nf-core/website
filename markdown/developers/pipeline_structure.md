# nf-core contribution guidelines

All nf-core pipelines adhere to common design principles that have been crafted through discussions and extensive real-world experiences with pipelines over time.
We understand, however, that not all about nf-core pipelines is obvious and can be confusing to people new to nf-core.
To help you along getting started with nf-core, writing your pipeline or contribute otherwise, we have created this documentation.
Hopefully you will find most things that might confuse you described here, and if not, be sure to bug the nf-core Slack channel!

## General

To start writing a new pipeline you should use the `nf-core create` command (see [docs](/developers/adding_pipelines)).
This will create a skeleton pipeline with a number of files, each with different purposes. Here's what they are for:

### Files

* `main.nf`: This is the main nextflow file which will get executed if the pipeline is run. It parses all parameters, loads the modules with options and defines the specific order in which these are executed.

* `nextflow.config`: The main nextflow configuration file. It contains the default pipeline parameters, nextflow configuration options and information like pipeline and minimum nextflow version, among others.
  The `nextflow.config` also defines different configuration profiles that can be used to run the pipeline. See the [Configuration docs](/usage/configuration) for more information.

* `README.md`: Basic information about the pipeline and usage

* `nextflow_json.schema`: The JSON schema file is used for pipeline parameter specification. This is automatically created using the `nf-core schema build` command. It is used for printing command-line help, validating input parameters, building the website docs and for building pipeline launch interfaces (web and cli).

* `CHANGELOG.md`: Information about the changes made to the pipeline for each release.

* `LICENSE`: The license - should be MIT

* `CODE_OF_CONDUCT.md`: The nf-core code of conduct.

* `.gitattributes`: Git settings, primarily getting the `.config` files to render with Nextflow syntax highlighting on <github.com>

* `.gitignore`: Files that should be ignored by git.

### Directories

* `.github/`: Other GitHub specific files, e.g. for specifying templates and GitHub actions

* `assets/`: Any additional files needed for the pipeline

* `bin/`: Directory for scripts that must be directly accessible within a pipeline process. Anything in this directory can be directly called from within Nextflow processes.

* `conf/`: Configuration files, including a `base.config` file which is always loaded into `nextflow.config` and describes basic pipeline configurations, like CPU and memory usage for processes with low, medium and high requirements. Additionaly, most pipelines also have a `igenomes.config` file which describes the locations of popular genomes that can be automatically downloaded for a pipeline run. Finally, the `test.config` and `test_full.config` files are test configurations that are loaded during test runs. Since DSL2, it also contains a `modules.config` file, which defines module-specific configurations and is explained further down in the "DSL2 and modules" section.

* `docs/`: Markdown files for documenting the pipeline

* `lib/`: The lib directory contains groovy utility functions. Currently, this includes the `Schema.groovy`, `Headers.groovy`, `Completion.groovy` and `Checks.groovy` files. These are called from within the nf-core pipeline to do common pipeline tasks, like formatting the header for the command line. Any larger groovy function you want to call in the code should therefore go in here.

* `modules/`: Contains pipeline-specific and common nf-core modules

## Continuous integration testing

To assure that nf-core pipelines don't break after some change is made to the code, we use automated continuous integration (CI) testing. This is done via GitHub actions, which are defined in the `.github/workflows` directory. Parameters and file paths are set in the `conf/test.config` and `conf/test_full.config`. Please see also [here](/developers/adding_pipelines#add-some-test-data) for how to set-up the test workflow for your pipeline.

## DSL2 and modules

Nextflow DSL2 allows for a more modularized design of pipelines and the reuse of components. Currently, most nf-core pipelines are still entirely written in DSL1, but in the near future all pipelines will be written in DSL2. The nf-core team has developed a set of design patterns on how to best implement DSL2 pipelines, which should be used by all nf-core pipelines in order to assure standardization and the reuse of components. The following is meant to help understand certain design choices and how a nf-core DSL2 pipeline should be build.

### Modules

Each pipeline has a `modules` directory which contains all the module code. A module here depicts a single process which involves - if possible - only a single tool/software. The `modules` directory is furthermore divided into `local`and `nf-core` sub-directories, which themselves each have a `process`/`software` and `subworkflow` directory. Modules contained in the `local` directory are specific to the pipeline, whereas `nf-core` modules are installed from the `nf-core/modules` repository. For instance, most pipelines that involve FastQ files will run the FastQC tool for quality control. The module needed for this can be easily reused from the `nf-core/modules` directory using the `nf-core/tools`package. The `process` directories contain modules which define single processes, which smaller workflows are contained in the `subworkflow` directories.

All modules load utility functions from a `functions.nf` script that must be contained in the `modules/local/process` directory. It contains simple functions to initialize the module options, get the software version, save files and get a path from a string. For further explanations of modules and how they should be structured in DSL2 pipelines, check out the [nf-core modules repo](https://github.com/nf-core/modules).

### Module parameters

One thing that might not be straightforward is how module parameters are handled in nf-core DSL2 pipelines. Every module and subworkflow, when loaded into the pipeline, has to be passed a groovy map containing module options. For single processes this is typically only a single `options` map, while for subworkflows these can be several maps that are then passed down to the correct processes within the subworkflows. These `options` maps are directly loaded from the `modules.config` file (contained in the pipeline `conf` directory), which is the place where all additional and optional parameters for modules are stored. Modules should be build in a way such that they are flexible with respect to the parameters, so that most command line parameters can be passed to them via the `modules.config`. This way, all command line parameters and other options can be modified within a single script, which makes it easy for users to adjust the pipeline and at the same time makes modules more reusable.

The `modules.config` file should contain a `params.modules` dictionary which lists every module used in the pipeline. For each module, the following fields can be specified:

* `args`: additional arguments appended to command in the module
* `args2`: Second set of arguments append to command in the module (multi-tool modules)
* `publish_dir`: Directory to publish the results
* `publish_by_id`: Publish results in separate folder by meta.id value
* `publish_files`: Groovy map where key = "file_ext" and value = "directory" to publish results for that file extension. The value of "directory" is appended to the standard "publish_dir" path as defined above. If publish_files == null (unspecified) all files are published. If publish_files == false no files are published.
* `suffix`: File name suffix for output files

### Sample meta information

In nf-core DSL2 pipelines, every channel that contains sample data in some way should also contain a `meta`variable, which must contain the fields `meta.id`, `meta.single_end` and `meta.strandedness`. The `meta` variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the `meta` variable. This meta information can easily be extracted from a samplesheet which specifies the input files. For an example process that reads a samplesheet, creates a `meta` variable and returns it along with the filepaths, have a look at the [rnaseq pipeline](https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf).
