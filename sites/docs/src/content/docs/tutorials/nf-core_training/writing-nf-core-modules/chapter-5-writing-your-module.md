---
title: "Chapter 5: Writing your module"
subtitle: "How to fill in the primary files of an nf-core module"
shortTitle: "Chapter 5: Writing"
---

## Introduction

Once you have generated the boilerplate template files for your module, you can start updating these to make your module function.
The boilerplate template files have lots of `TODO` comments that can help you give pointers how to write the module, as well example contents.

In this chapter we will go through each one, explaining what each section of each file is doing and why some bits are the way they are, as defined by the nf-core specifications.

## The `environment.yml` file

The first file that is generated is an conda `environment.yml` file.

If the `nf-core modules create` command was successfully able to find the tool within Bioconda, it should already have specified the channel, tool name, and the specific version of the tool, as you can see in the example below.

```yaml
---
# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
channels:
  - conda-forge
  - bioconda
dependencies:
  - "bioconda::drep=3.5.0"
```

You very rarely need to modify or touch this file.

The tool and version specified in this file will in most cases have a corresponding biocontainer (specified in the `main.nf` file).

You only need to modify the `environment.yaml` if your tool does not use tools from Bioconda (e.g. a more generic tool that is only on [`conda-forge`](https://conda-forge.org/)), and there is not an existing Bioconda recipe for the tool.

However, the nf-core guidelines strongly recommend conda support, as conda is the most accessible software management system.
So if your tool is not on Bioconda or conda-forge, we strongly recommend you add it to the appropriate repository.

:::info{title="Behind the scenes" collapse}
nf-core uses a separate conda file rather than defining within the `main.nf` script file's conda directive, for two main reasons:

- It will facilitate more [automated container building with Seqera's `wave` infrastructure](https://nf-co.re/blog/2024/seqera-containers-part-1)
- It makes multi-tool environments much easier to read and manage (vs. having multiple conda declarations on a single line)

  :::

## The `main.nf` file

This is the main file of the nf-core module, where the Nextflow module code itself is defined.

When you first generate the file, some parts of the module will already be filled in for you, and there will be many `TODO` comments to help guide you through writing the module.

There are 9 main Nextflow process blocks of an nf-core module

- `tag:`
- `label:`
- `conda:`
- `container:`
- `input:`
- `output:`
- `when:`
- `script:`
- `stub:`

:::info{title="Click here to see full 'raw' file example" collapse}

The boilerplate TODO comments have been removed for readability.

```nextflow
process DREP_COMPARE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0':
        'biocontainers/drep:3.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
```

:::

We will now go through each block of the module, and describe what each section is doing and summarise the general specification of each section.

### The `tag` block

The `tag` block is automatically propagated with a standard bit of metadata that is associated with each input channel containing a [meta map](https://nf-co.re/docs/contributing/components/meta_map).

```nextflow
tag "$meta.id"
```

If it is included, you generally do not need to modify this, as it assumed all nf-core pipelines will have the `id` element in each meta map and this can be changed by a pipeline developer in the pipeline code.

### The `label` block

The `label` block is automatically propagated by the option you specified during the boilerplate file generation.

```nextflow
label 'process_single'
```

nf-core defines a preset set of labels with associated default resources, as specified [here](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L22-L54).
You should have selected from one of these standard labels that approximately matches the typical requirements of the tool you are creating the module for.
If you change your mind after the initial boilerplate file generation, you can always modify it to one of the other standard labels.

You should not use custom label names to specify other 'custom' requirements if you feel your module's typical resource requirements do not match those of set in the nf-core pipeline template's standard labels.
The default resource specifications that are set for the standard labels can be modified by a pipeline developer when developing a pipeline level (in nf-core pipelines, this happens in the `conf/base.config` file).
Thus, you should stick with the standard labels, and let pipeline developers 'tune' the modules resources within the pipeline.
Furthermore a user can also override and refine these resource requirements even further for specific modules within a custom custom config file, by specifying in a `process` block by the module name with `withName:`.

This flexibility at the pipeline level means that it is unnecessary to use any non-standard labels, and provides the benefit of 'consistency' between all nf-core modules and nf-core pipelines on how each standard resource is defined.

### The `conda` block

The `conda` block tells Nextflow to use the associated `environment.yml` file to generate the required conda environment for this module, when a user has specified to use conda for running a pipeline.

```nextflow
conda "${moduleDir}/environment.yml"
```

You should not modify this, except in the very rare exceptional cases where a given tool does not support conda (typically proprietary tools), in which case you may remove the line.

### The `container` block

The `container` block tells Nextflow where to pull the tool's docker and singularity containers from.

```nextflow
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0':
    'biocontainers/drep:3.5.0--pyhdfd78af_0' }"
```

In most cases the `nf-core modules create` command should have auto-detected the name and versions of each container based on the tool's name, version, and build of the given container on Bioconda/biocontainers, as with the conda `environment.yml` file.

You should not modify this, except in the very rare exceptional cases where a given tool does not have a Bioconda recipe nor biocontainer (typically proprietary tools), in which case you may manually modify these lines to point to the corresponding container.

### The `input` block

The `input` block is where you will define what inputs the module will receive.

```nextflow
input:
tuple val(meta), path(bam)
```

The boilerplate template comes with an example of a common bioinformatic file format (`.bam`) and an associated [`meta` map](https://nf-co.re/docs/contributing/components/meta_map) (the meta map may be missing if you did not specify this during the boilerplate template file creation).

You should edit this section of the `main.nf` to define all the files your module will require.

The input block of an nf-core module should include **all possible files** via a `path()` declaration for each input file a given tool or tool subcommand can accept - this covers both mandatory **and** optional files.
You should make sure to check the documentation of your tool to ensure you capture all of these!

:::info{title="Behind the scenes" collapse}

We require _all_ files to be specified here because Nextflow processes _must_ receive files via this route to correctly 'stage' the files in a jobs working directory.

Furthermore, you cannot assume that your use case of the module will cover what other pipeline developers may need to use it for, thus you should try to be as comprehensive as possible.

While there is no 'native' Nextflow way to specify optional inputs, a pipeline developer can instead pass an empty list `[]` to the module to indicate that input is not being provided.
:::

In addition to files, mandatory tool argument/parameters (non-file) should also be specified here as a `val` input channel (in most cases without a meta map).

```nextflow
input:
tuple val(meta), path(bam)
val mode
```

These are often required arguments without a default, such as a `--mode` parameter.
You do not need to define dedicated channels for optional parameters or arguments, as these are specified via different mechanism described later.

:::info{title="Behind the scenes" collapse}
We also require all mandatory non-file arguments to ensure that a particular module has all the _minimum_ required information for a tool to execute, so it can run 'out of the box' without any extra configuration (see `ext.args` below).
:::

Input channels names generally follow the convention of using the file format as the variable name or the tool's argument name that is used in the command to refer to the input object.

If the tool has multiple input files, you should make one channel per input file.
Each new channel can have its own meta map, with each subsequent one being called `meta2`, `meta3`.
Generally we recommend all input channels should have a `meta`, except if you can guarantee that you would never have more than a single version of given input file (e.g., a configuration file that applies to every single sample/input data file).

In a few exceptional cases, multiple very closely related files can be specified 'inline' with a single tuple channel, such as index files with the main file (e.g. `.bai` file with a `.bam` file), where the index file cannot be used by itself without the main file.

```nextflow
input:
tuple val(meta), path(bam), path(bai)
```

:::info{title="Behind the scenes" collapse}
The purposes of the 'one channel = one file' specification is to ensure as much consistency as possible for pipeline developers.

If we allowed high heterogeneity in the number and ordering of files in a single input channel, it makes it harder for pipeline developers to know how to prepare input/output channels.

For a pipeline developer, they can use Nextflow's [`.multiMap()`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) operator to split their own multi-element channel into a separate 'sub-channels' to ensure the objects in the elements stay associated with each other even though being separated into channels.
:::

:::tip{title="Examples" collapse}

An example of an `input:` block of a tool that requires a single input file and a mandatory argument:

```nextflow
input:
tuple val(meta), path(fasta)
val(output_format)
```

An example of an input block that requires multiple input files:

```nextflow
input:
tuple val(meta), path(reads)
tuple val(meta2), path(fasta)
tuple val(meta3), path(fasta_index)
```

An example of an input block that requires multiple input files, where the first channel's file has a closely associated index file:

```nextflow
input:
tuple val(meta), path(input), path(intervals)
path  fasta
```

:::

### The `output` block

The `output:` block is where you will define all possible output files a tool can produce.

```nextflow
output:
path "versions.yml"           , emit: versions
```

The boilerplate template code by default includes a mandatory `versions.yml` file emission.
This is because all nf-core pipelines must report the versions of every tool used within the pipeline.

Above the `versions.yml` output channel line you can include the output definitions of all files your tool(/subcommand) will produce.
You should already define as many output files as possible (whether produced by default, or optionally), to make as much functionality available to future users as possible.
If your module uses a meta map, you should emit all files with a meta map, except for `versions.yml`.

```nextflow
output:
tuple val(meta), path("${prefix}.bam"),  emit: bam,  optional: true
tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
tuple val(meta), path("${prefix}.sam"),  emit: sam,  optional: true
path  "versions.yml",
```

Each output file type should have it's own output `emit` entry, e.g. one per file format.
In some cases if two mutually exclusive alternative file formats can be created by the tool, but serve the same purpose (such as two different indexing formats, for example `.bai` and `.csi` for `.bam files`), these can be captured in the same pattern and output channel.

:::info{title="Behind the scenes" collapse}
We often will recommend including as many input/output files as possible in the first creation of a tool.
This is because want to both make sure your nf-core module can be used by many people as possible, but also to minimise the number of channel structure changes/updates a module has to go through.

By being as comprehensive as possible during initial creation, we maintain stability and familiarity of how the module is used by pipeline developers and thus reduces the risk of breaking pipelines when they update modules.
:::

Generally, we recommend that where possible all output files should be compressed (if a tool can accept compressed inputs), thus you should ensure the path patterns take this into account.

```nextflow
output:
tuple val(meta), path("*.{vcf,vcf.gz}"),  emit: vcf
```

The `emit:` channel names of each output file should generally correspond to the file format/suffix of each format.
If an output file is not always produced, this should be indicated as such using Nextflow's `optional: true` option.

:::tip{title="Examples" collapse}
An example of a module that produces multiple files, some of which are optional

```bash
output:
tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
tuple val(meta), path('*.json')           , emit: json
tuple val(meta), path('*.html')           , emit: html
tuple val(meta), path('*.log')            , emit: log
tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
path "versions.yml"                       , emit: versions
```

An example of a module that produces three mutually exclusive primary output files:

```nextflow
output:
tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
tuple val(meta), path("*.csi")                    , emit: csi, optional: true
path "versions.yml"                               , emit: versions
```

An example of a module that produces multiple files all of which can support compressed variants:

```nextflow
output:
tuple val(meta), path('*.{blast,blast.gz}'), optional: true, emit: blast
tuple val(meta), path('*.{xml,xml.gz}')    , optional: true, emit: xml
tuple val(meta), path('*.{txt,txt.gz}')    , optional: true, emit: txt
tuple val(meta), path('*.{daa,daa.gz}')    , optional: true, emit: daa
tuple val(meta), path('*.{sam,sam.gz}')    , optional: true, emit: sam
tuple val(meta), path('*.{tsv,tsv.gz}')    , optional: true, emit: tsv
tuple val(meta), path('*.{paf,paf.gz}')    , optional: true, emit: paf
path "versions.yml"                        , emit: versions
```

:::

### The `when` block

The `when` block can be used by pipeline developers to dynamically activate or deactivate the execution of the given module.

```nextflow
when:
task.ext.when == null || task.ext.when
```

You should not modify or remove this block.

### The `script` block

The `script` block is the main part of the module, where the actual command that the module will execute within pipelines is defined.
You will do most modification of your `main.nf` here.

```nextflow
script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
samtools \\
    sort \\
    $args \\
    -@ $task.cpus \\
    -o ${prefix}.bam \\
    -T $prefix \\
    $bam

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    drep: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
"""
```

The boilerplate template already comes with an example samtools command for you, as well a second versions command creating the `versions.yml` HEREDOC, both of which you will need to replace for your respective tool.

At the top of the block you have two standard variables that should not be removed.

- 📛 The`args` variable is how a pipeline developer can 'inject' optional parameters into the command itself.
  - The `ext.args` comes from the process scope of a nextflow configuration file - in nf-core pipelines this is defined in the `modules.config`
  - More information on this is described in the 'using in pipelines' chapter
- 📛 The `prefix` variable is used for the default output file basename that all output files should use
  - The default corresponds to the default `id` value of the primary input channel's meta map that all nf-core pipelines use.

Within the script section, you can replace the example with your own command, split across multiple lines with escaped backslashes for readability (see example above), specifying all the variables you need based on your input channels.

All commands must have the `$args` variable used within it, and if supported, the number of CPUs the command can use.

You can add additional variables after the `args` and `prefix` if you need to dynamically inject input files in different ways (e.g., in genomics if you have to use different flags for specifying single- or paired-end sequencing data files)

```nextflow
script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def single_end  = meta.single_end ? "--single" : ""
"""
elprep merge \\
       input/ \\
       output/${prefix}.bam \\
       $args \\
       ${single_end} \\
...
"""
```

You must keep the `versions.yml` HEREDOC, but replace the command to generate a clean 'version' (i.e. just `1.2.3`, without leading `v` etc) of the version from the tools `--version` output.
We typically use standard UNIX command line tools for such a clean up, such as `cut` or `sed` for this purpose.
If your tool does not have way of emitting this information, you can use specify a dedicated variable after the `args` and `prefix` variables, that can then be specified and interpreted within the HEREDOC.

```nextflow
script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def VERSION='2.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping
"""
...

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    scimap: $VERSION
END_VERSIONS
```

While there is a relatively strict 'one module == one command/tool' rule, in some exceptional cases you may need to chain multiple tools (e.g. piping) or have a second command (e.g. compression).
In these cases you must define a second `args` variable with `args2`, `args3` and so on for each part of the pipe or tool.

It is important that you do not use any custom `meta` elements in modules!

:::info{title="Behind the scenes" collapse}
Custom meta fields will not be standardised across tools, and add extra burden on pipeline developers who may want to encode the information in their meta maps by different names.

All optional/additional information for dynamic construction of a command away from the bare minimum execution should be left to the pipeline developer to inject via the `ext.args`
:::

:::warning
We do not provide any explicit examples here as there is a lot of variation across different tools.
We recommend you to just check multiple modules within the [nf-core/modules GitHub repository](https://github.com/nf-core/modules) to get further guidance.
:::

### The `stub` block

The `stub` block is meant simulate the output emissions of a module during a `-dry-run` of a pipeline.

```nextflow
  stub:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  touch ${prefix}.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      drep: \$(samtools --version |& sed '1!d ; s/samtools //')
  END_VERSIONS
  """
```

For nf-core modules, the stub block should essentially copy the `script` block.
However instead of executing the actual command, you just create empty files with the same names as actual output files the tool would create using one or more `touch` commands.
If your tool outputs gzipped files, you should pipe the output of an empty echo command into `gzip`,

```nextflow
stub:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
echo "" | gzip > ${prefix}.txt.gz

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
"""
```

You should create files that can be picked up by every `output` channel.

If you are editing in VScode, and get [Nextflow language server](https://www.nextflow.io/docs/latest/vscode.html) errors, you can either copy the command from the `script` block and wrap it in an `echo`, or you may simplify remove the `args` definition.

```nextflow
stub:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
echo "$args"
echo "" | gzip > ${prefix}.txt.gz

cat <<-END_VERSIONS > versions.yml
"${task.process}":
tool: \$(samtools --version |& sed '1!d ; s/samtools //')
END_VERSIONS
"""
```

It is important to still copy/update the `versions.yml` HEREDOC command from the `script:` block, as this will still be executed during the dry run.

## The `meta.yaml` file

Another 'unique' component of an nf-core module is the `meta.yaml` file.

The purpose of this file is to better document all parts of the module, with additional descriptions, keywords, and links to the tool's original resources.
Much of this information is used to help make the modules searchable on the nf-core website [modules page](https://nf-co.re/modules), make the wrapper findable on external websites and databases, and to prepare the community's move towards automated linkage between modules.

The `meta.yaml` can be split into four main sections

- Tool and module description
- Input channel descriptions
- Output channel descriptions
- Contributors

:::info{title="Click here to see full 'raw' file example" collapse}

The boilerplate TODO comments have been removed for readability.

```yaml
---
# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/meta-schema.json
name: "drep_compare"
description: write your description here
keywords:
  - sort
  - example
  - genomics
tools:
  - "drep":
      description: "De-replication of microbial genomes assembled from multiple samples"
      homepage: "https://drep.readthedocs.io/en/latest/"
      documentation: "https://drep.readthedocs.io/en/latest/"
      tool_dev_url: "https://github.com/MrOlm/drep"
      doi: ""
      licence: ["MIT"]
      identifier: biotools:drep

input:
  - - meta:
        type: map
        description: |
          Groovy Map containing sample information
          e.g. `[ id:'sample1', single_end:false ]`
    - bam:
        type: file
        description: Sorted BAM/CRAM/SAM file
        pattern: "*.{bam,cram,sam}"
        ontologies:
          - edam: "http://edamontology.org/format_25722"
          - edam: "http://edamontology.org/format_2573"
          - edam: "http://edamontology.org/format_3462"

output:
  - bam:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. `[ id:'sample1', single_end:false ]`
      - "*.bam":
          type: file
          description: Sorted BAM/CRAM/SAM file
          pattern: "*.{bam,cram,sam}"
          ontologies:
            - edam: "http://edamontology.org/format_25722"
            - edam: "http://edamontology.org/format_2573"
            - edam: "http://edamontology.org/format_3462"

  - versions:
      - "versions.yml":
          type: file
          description: File containing software versions
          pattern: "versions.yml"

authors:
  - "@jfy133"
maintainers:
  - "@jfy133"
```

:::

We will now go through each of these sections.

### Tool and module description

This section provides metadata about the tool and the module itself to make it more findable e.g. on the nf-core website.

```yaml
name: "drep_compare"
description: write your description here
keywords:
  - sort
  - example
  - genomics
tools:
  - "drep":
      description: "De-replication of microbial genomes assembled from multiple samples"
      homepage: "https://drep.readthedocs.io/en/latest/"
      documentation: "https://drep.readthedocs.io/en/latest/"
      tool_dev_url: "https://github.com/MrOlm/drep"
      doi: ""
      licence: ["MIT"]
      identifier: biotools:drep
```

Some parts of this section are automatically filled in for you based on the name of the module, and information already recorded in metadata of the Bioconda recipe that is used in the previously described `environment.yaml`

In most cases, you only need to update the keywords, and the module description.

For the module description (i.e., second line after the first `name:` line), you need to describe the purpose of the _module_ in a short sentence.

For the keywords, you need to have a minimum of 3 - all lower case - and should include things such as the tool's name, subcommand, words to describe the 'action' the module is performing such as 'sort' or 'filter', but can also be words describing the research field the module used - such as genomics, metagenomics, earth sciences, and the file formats the action is used on.
Ultimately, the more keywords the better.

Sometimes you may need to fill in the URLs to the tool's source code and documentation.
It is OK to leave these empty or duplicate URLs if the tool does not have a dedicated documentation page, for example.

Note that if you are using more than one tool in your module (e.g. the tool and gzip), you will need to fill in the tool description URLs for each tool.

### Input channel descriptions

This section describes each input channel of the module.

```yaml
input:
  # Only when we have meta
  - - meta:
        type: map
        description: |
          Groovy Map containing sample information
          e.g. `[ id:'sample1', single_end:false ]`

    - bam:
        type: file
        description: Sorted BAM/CRAM/SAM file
        pattern: "*.{bam,cram,sam}"
        ontologies:
          - edam: "http://edamontology.org/format_25722"
          - edam: "http://edamontology.org/format_2573"
          - edam: "http://edamontology.org/format_3462"
```

In the default example added by the boilerplate template, you see a single input channel described - one with a meta map and a single path element.

In this section you need to update the the input channel name, the description, pattern, and ontologies.

The type follows a fixed set of [categories](https://nf-co.re/docs/guidelines/components/modules#input-and-output-channel-types), such as `file`, `integer`, `boolean`.

The description should be descriptive in what they contain or how they should be prepared, e.g. 'A TSV file containing 5 columns generated by XYZ', not just a `TSV`.

The pattern should match the corresponding input expected by the tool itself.

The ontologies section should provide links to respective entry matching the input type within the the controlled vocabulary of the (typically) [EDAM ontology](https://www.ebi.ac.uk/ols4/ontologies/edam).
It is up to you how specific to be - whether extremely specific (fastq), or generic (text file).

Note that you will have to make an entry for each of your input channels.
If you copy and paste, make sure to edit each section to be unique for each of your input channels.
When each of the additional meta files have their own `meta` maps, the name should be updated to `meta2` etc accordingly.

:::tip
You can use the `nf-core modules lint` command (described later) to automatically pick up all the input channels and put placeholders elements in for you
:::

:::tip{title="Examples" collapse}
Example of a complete `meta.yaml` using ontologies and multiple meta entries (see below for `outputs:`)

```yml
name: bwamem2_mem
description: Performs fastq alignment to a fasta reference using BWA
keywords:
  - mem
  - bwa
  - alignment
  - map
  - fastq
  - bam
  - sam
tools:
  - bwa:
      description: |
        BWA-mem2 is a software package for mapping DNA sequences against
        a large reference genome, such as the human genome.
      homepage: https://github.com/bwa-mem2/bwa-mem2
      documentation: http://www.htslib.org/doc/samtools.html
      arxiv: arXiv:1303.3997
      licence: ["MIT"]
      identifier: "biotools:bwa-mem2"
input:
  - - meta:
        type: map
        description: |
          Groovy Map containing sample information
          e.g. [ id:'test', single_end:false ]
    - reads:
        type: file
        description: |
          List of input FastQ files of size 1 and 2 for single-end and paired-end data,
          respectively.
        ontologies:
          - edam: "http://edamontology.org/data_2044" # Sequence
          - edam: "http://edamontology.org/format_1930" # FASTQ
  - - meta2:
        type: map
        description: |
          Groovy Map containing reference/index information
          e.g. [ id:'test' ]
    - index:
        type: file
        description: BWA genome index files
        pattern: "Directory containing BWA index *.{0132,amb,ann,bwt.2bit.64,pac}"
        ontologies:
          - edam: "http://edamontology.org/data_3210" # Genome index
  - - meta3:
        type: map
        description: |
          Groovy Map containing reference information
          e.g. [ id:'genome' ]
    - fasta:
        type: file
        description: Reference genome in FASTA format
        pattern: "*.{fa,fasta,fna}"
        ontologies:
          - edam: "http://edamontology.org/data_2044" # Sequence
          - edam: "http://edamontology.org/format_1929" # FASTA
  - - sort_bam:
        type: boolean
        description: use samtools sort (true) or samtools view (false)
        pattern: "true or false"
output:
  - sam:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. [ id:'test', single_end:false ]
      - "*.sam":
          type: file
          description: Output SAM file containing read alignments
          pattern: "*.{sam}"
          ontologies:
            - edam: "http://edamontology.org/format_2573" # SAM
  - bam:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. [ id:'test', single_end:false ]
      - "*.bam":
          type: file
          description: Output BAM file containing read alignments
          pattern: "*.{bam}"
          ontologies:
            - edam: "http://edamontology.org/format_2572" # BAM
  - cram:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. [ id:'test', single_end:false ]
      - "*.cram":
          type: file
          description: Output CRAM file containing read alignments
          pattern: "*.{cram}"
          ontologies:
            - edam: "http://edamontology.org/format_3462" # CRAM
  - crai:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. [ id:'test', single_end:false ]
      - "*.crai":
          type: file
          description: Index file for CRAM file
          pattern: "*.{crai}"
  - csi:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. [ id:'test', single_end:false ]
      - "*.csi":
          type: file
          description: Index file for BAM file
          pattern: "*.{csi}"
  - versions:
      - versions.yml:
          type: file
          description: File containing software versions
          pattern: "versions.yml"
authors:
  - "@maxulysse"
  - "@matthdsm"
maintainers:
  - "@maxulysse"
  - "@matthdsm"
```

:::

### Output channel descriptions

This section follows the same concept as the input section in the section above.

```yaml
output:
  - bam:
      #Only when we have meta
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. `[ id:'sample1', single_end:false ]`
      - "*.bam":
          type: file
          description: Sorted BAM/CRAM/SAM file
          pattern: "*.{bam,cram,sam}"
          ontologies:
            - edam: "http://edamontology.org/format_25722"
            - edam: "http://edamontology.org/format_2573"
            - edam: "http://edamontology.org/format_3462"
```

The only difference is that the 'key' of output _file_ (not the channel) should match that of the pattern used in the output channel as in the `main.nf`.

:::tip
You can use the `nf-core modules lint` command (described later) to automatically pick up all the output channels and put placeholders elements in for you
:::

As with the the `input:` block you need to update the type, description, pattern, and ontologies.

### Contributors

This section simply describes the original author and other contributors who have added or updated the module since original creation.

```yaml
authors:
  - "@jfy133"
maintainers:
  - "@jfy133"
```

You do not have to modify this section as it is automatically propagated during the boilerplate template file creation.

If you later update another module originally by someone else, you are welcome to add yourself to the 'maintainers:' section.

## Summary

In this chapter we have gone through each file in fine detail, describing how each part of each file should be filled out.
This tutorial cannot cover every single point of the specifications, nor cover every edge case.
We therefore highly recommend to always refer to the [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules) when writing a module.

In the next chapter we will describe the other half of an nf-core modules - the files used for unit testing of the modules.
