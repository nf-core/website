---
title: DSL2 Modules
subtitle: Documentation and links for nf-core modules
---

If you decide to upload a module to `nf-core/modules` then this will ensure that it will become available to all nf-core pipelines, and to everyone within the Nextflow community! See [`modules/`](https://github.com/nf-core/modules/tree/master/modules) for examples.

See [the dsl2 modules tutorial](https://nf-co.re/developers/tutorials/dsl2_modules_tutorial) for how to add a module!

## Terminology

The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components.

### Module

A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. At present, this repository has been created to only host atomic module files that should be added to the [`modules/`](https://github.com/nf-core/modules/tree/master/modules) directory of nf-core/modules along with the required documentation and tests.

### Sub-workflow

A chain of multiple modules that offer a higher-level of functionality within the context of a pipeline. For example, a sub-workflow to run multiple QC tools with FastQ files as input. Sub-workflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. As it stands, this repository will not host sub-workflows although this may change in the future since well-written sub-workflows will be the most powerful aspect of DSL2.

### Workflow

What DSL1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL1, or by using a combination of DSL2 individual modules and sub-workflows.

## Guidelines

### General

- All non-mandatory command-line tool non-file arguments MUST be provided as a string via the `$args` variable, which is assigned to using the `task.ext.args` variable. The value of `task.ext.args` is supplied from the `modules.config` file by assigning a string value to `ext.args`.

    `<module>.nf`:

    ```nextflow
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc \\
        $args \\
         <...>
    """
    ```

    `modules.config`:

    ```nextflow
    process {
        withName: <module> {
            ext.args = [                                                          // Assign either a string, closure which returns a string
                '--quiet',
                params.fastqc_kmer_size ? "-k ${params.fastqc_kmer_size}" : ''    // Parameter dependent values can be provided like so
            ].join(' ')                                                           // Join converts the list here to a string.
            ext.prefix = { "${meta.id}" }                                         // A closure can be used to access variables defined in the script
        }
    }
    ```

- Software that can be piped together SHOULD be added to separate module files
unless there is a run-time, storage advantage in implementing in this way. For example,
using a combination of `bwa` and `samtools` to output a BAM file instead of a SAM file:

    ```bash
    bwa mem | samtools view -B -T ref.fasta
    ```

- Where applicable, the usage and generation of compressed files SHOULD be enforced as input and output, respectively:
  - `*.fastq.gz` and NOT `*.fastq`
  - `*.bam` and NOT `*.sam`

  If a tool does not support compressed input or output natively, we RECOMMEND passing the
  uncompressed data via unix pipes, such that it never gets written to disk, e.g.

  ```bash
  gzip -cdf $input | tool | gzip > $output
  ```

  The `-f` option makes `gzip` auto-detect if the input is compressed or not.

  If a tool cannot read from STDIN, or has multiple input files, it is possible to use
  named pipes:

  ```bash
  mkfifo input1_uncompressed input2_uncompressed
  gzip -cdf $input1 > input1_uncompressed &
  gzip -cdf $input2 > input2_uncompressed &
  tool input1_uncompressed input2_uncompressed > $output
  ```

  Only if a tool reads the input multiple times, it is required to uncompress the
  file before running the tool.

- Where applicable, each module command MUST emit a file `versions.yml` containing the version number for each tool executed by the module, e.g.

    ```bash
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
        samtools: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$// )
    END_VERSION
    ```

    resulting in, for instance,

    ```yaml
    "FASTQC":
        fastqc: 0.11.9
        samtools: 1.12
    ```

    All reported versions MUST be without a leading `v` or similar, i.e. must start with a numeric character.

    We chose a [HEREDOC](https://tldp.org/LDP/abs/html/here-docs.html) over piping into the versions file
    line-by-line as we believe the latter makes it easy to accidentally overwrite the file. Moreover, the exit status
    of the sub-shells evaluated in within the HEREDOC is ignored, ensuring that a tool's version command does
    not erroneously terminate the module.

    If the software is unable to output a version number on the command-line then a variable called `VERSION` can be manually specified to create this file e.g. [homer/annotatepeaks module](https://github.com/nf-core/modules/blob/master/modules/homer/annotatepeaks/main.nf). If the HEREDOC cannot be used because the script is not bash, the versions.yml must be written directly e.g. [ascat module](https://github.com/nf-core/modules/blob/master/modules/ascat/main.nf).

- The process definition MUST NOT change the `when` statement. `when` conditions can instead be supplied using the `process.ext.when` directive in a configuration file.

    ```groovy
    process {
        withName: 'FOO' {
            ext.when = !params.skip_module
        }
        withName: 'BAR' {
            ext.when = { meta.single_end }
        }
    }
    ```

### Naming conventions

- The directory structure for the module name must be all lowercase e.g. [`modules/bwa/mem/`](https://github.com/nf-core/modules/tree/master/modules/bwa/mem/). The name of the software (i.e. `bwa`) and tool (i.e. `mem`) MUST be all one word.

- The process name in the module file MUST be all uppercase e.g. `process BWA_MEM {`. The name of the software (i.e. `BWA`) and tool (i.e. `MEM`) MUST be all one word separated by an underscore.

- All parameter names MUST follow the `snake_case` convention.

- All function names MUST follow the `camelCase` convention.

- Output file (and/or directory) names SHOULD just consist of only `${prefix}` and the file-format suffix (e.g. `${prefix}.fq.gz` or `${prefix}.bam`).
  - This is primarily for re-usability so that other developers have complete flexibility to name their output files however they wish when using the same module.
  - As a result of using this syntax, if the module has the same named inputs and outputs then you can add a line in the `script` section like below (another example [here](https://github.com/nf-core/modules/blob/e20e57f90b6787ac9a010a980cf6ea98bd990046/modules/lima/main.nf#L37)) which will raise an error asking the developer to change the `args.prefix` variable to rename the output files so they don't clash.
  
    ```nextflow
    script:
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    ```

### Input/output options

- Input channel declarations MUST be defined for all _possible_ input files (i.e. both required and optional files).
  - Directly associated auxiliary files to an input file MAY be defined within the same input channel alongside the main input channel  (e.g. [BAM and BAI](https://github.com/nf-core/modules/blob/e937c7950af70930d1f34bb961403d9d2aa81c7d/modules/samtools/flagstat/main.nf#L22)).
  - Other generic auxiliary files used across different input files (e.g. common reference sequences) MAY be defined using a dedicated input channel (e.g. [reference files](https://github.com/nf-core/modules/blob/3cabc95d0ed8a5a4e07b8f9b1d1f7ff9a70f61e1/modules/bwa/mem/main.nf#L21-L23)).

- Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

- Optional inputs are not currently supported by Nextflow. However, passing an empty list (`[]`) instead of a file as a module parameter can be used to work around this issue.

### Module parameters

- A module file SHOULD only define input and output files as command-line parameters to be executed within the process.

- All `params` within the module MUST be initialised and used in the local context of the module. In other words, named `params` defined in the parent workflow MUST NOT be assumed to be passed to the module to allow developers to call their parameters whatever they want. In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

- If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

- Any parameters that need to be evaluated in the context of a particular sample e.g. single-end/paired-end data MUST also be defined within the process.

### Resource requirements

- An appropriate resource `label` MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config#L29-L46) e.g. `process_low`, `process_medium` or `process_high`.

- If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

- If a module contains _multiple_ tools that supports multi-threading (e.g. [piping output into a samtools command](https://github.com/nf-core/modules/blob/28b023e6f4d0d2745406d9dc6e38006882804e67/modules/bowtie2/align/main.nf#L32-L46)), you MUST assign cpus per tool such that the total number of used CPUs does not exceed `task.cpus`.
  - For example, combining two (or more) tools that both (all) have multi-threading, this can be assigned to the variable [`split_cpus`](https://github.com/nf-core/modules/blob/28b023e6f4d0d2745406d9dc6e38006882804e67/modules/bowtie2/align/main.nf#L32)
  - If one tool is multi-threaded and another uses a single thread, you can specify directly in the command itself e.g. with [`${task.cpus - 1}`](https://github.com/nf-core/modules/blob/6e68c1af9a514bb056c0513ebba6764efd6750fc/modules/bwa/sampe/main.nf#L42-L43)

### Software requirements

[BioContainers](https://biocontainers.pro/#/) is a registry of Docker and Singularity containers automatically created from all of the software packages on [Bioconda](https://bioconda.github.io/). Where possible we will use BioContainers to fetch pre-built software containers and Bioconda to install software using Conda.

- Software requirements SHOULD be declared within the module file using the Nextflow `container` directive. For single-tool BioContainers, the `nf-core modules create` command will automatically fetch and fill-in the appropriate Conda / Docker / Singularity definitions by parsing the information provided in the first part of the module name:

    ```nextflow
    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"
    ```

- If the software is available on Conda it MUST also be defined using the Nextflow `conda` directive. Using `bioconda::bwa=0.7.17` as an example, software MUST be pinned to the channel (i.e. `bioconda`) and version (i.e. `0.7.17`). Conda packages MUST not be pinned to a build because they can vary on different platforms.

- If required, multi-tool containers may also be available on BioContainers e.g. [`bwa` and `samtools`](https://biocontainers.pro/#/tools/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40). You can install and use the [`galaxy-tool-util`](https://anaconda.org/bioconda/galaxy-tool-util) package to search for both single- and multi-tool containers available in Conda, Docker and Singularity format. e.g. to search for Docker (hosted on Quay.io) and Singularity multi-tool containers with both `bowtie` and `samtools` installed you can use the following command:

    ```console
    mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep "mulled"
    ```

    > NB: Build information for all tools within a multi-tool container can be obtained in the `/usr/local/conda-meta/history` file within the container.

- It is also possible for a new multi-tool container to be built and added to BioContainers by submitting a pull request on their [`multi-package-containers`](https://github.com/BioContainers/multi-package-containers) repository.
  - Fork the [multi-package-containers repository](https://github.com/BioContainers/multi-package-containers)
  - Make a change to the `hash.tsv` file in the `combinations` directory see [here](https://github.com/aunderwo/multi-package-containers/blob/master/combinations/hash.tsv#L124) for an example where `pysam=0.16.0.1,biopython=1.78` was added.
  - Commit the code and then make a pull request to the original repo, for [example](https://github.com/BioContainers/multi-package-containers/pull/1661)
  - Once the PR has been accepted a container will get built and you can find it using  a search tool in the `galaxy-tool-util conda` package

      ```console
      mulled-search --destination quay singularity conda  --search pysam biopython  | grep "mulled"
      quay         mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  docker pull quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
      singularity  mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  wget https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
      ```

    - You can copy and paste the `mulled-*` path into the relevant Docker and Singularity lines in the Nextflow `process` definition of your module
    - To confirm that this is correct. Spin up a temporary Docker container

      ```console
      docker run --rm -it quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0  /bin/sh
      ```

      And in the command prompt type

      ```console
      $ grep specs /usr/local/conda-meta/history
      # update specs: ['biopython=1.78', 'pysam=0.16.0.1']
      ```

      The packages should reflect those added to the multi-package-containers repo `hash.tsv` file

- If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. We will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).

### Publishing results

Fomerly, results were published using a custom `publishDir` definition, customised using a Groovy Map defined by `params.modules`. This system has been replaced using Nextflow's native [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) defined directly in a pipeline workflow's `modules.config` (see [here](https://github.com/nf-core/rnaseq/blob/f7702d5b76a1351e2e7796a5ed3f59943a139fbf/conf/modules.config#L100-L106) for a simple example)

### Test data config file

If a new test dataset is added to [`tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), check that the config name of the added file(s) follows the scheme of the entire file name with dots replaced with underscores.

For example: the nf-core/test-datasets file `genomics/sarscov2/genome/genome.fasta` labelled as `genome_fasta`, or `genomics/sarscov2/genome/genome.fasta.fai` as `genome_fasta_fai`.

### Using a stub test when required test data is too big

If the module absolute cannot run using tiny test data, there is a possibility to add [stub-run](https://www.nextflow.io/docs/edge/process.html#stub) to the test.yml. In this case it is required to test the module using larger scale data and document how this is done. In addition, an extra script-block labeled `stub:` must be added, and this block must create dummy versions of all expected output files as well as the `versions.yml`. An example is found in the [ascat module](https://github.com/nf-core/modules/blob/master/modules/ascat/main.nf). In the `test.yml` the `-stub-run` argument is written as well as the md5sums for each of the files that are added in the stub-block. This causes the stub-code block to be activated when the unit test is run ([example](https://github.com/nf-core/modules/blob/master/tests/modules/ascat/test.yml)):

```console
nextflow run tests/modules/<nameofmodule> -entry test_<nameofmodule> -c tests/config/nextflow.config -stub-run
```

## Help

For further information or help, don't hesitate to get in touch on [Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)).
