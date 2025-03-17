---
title: "Chapter 5: Writing your module"
subtitle: "How to fill in the primary files of an nf-core module"
shortTitle: "Chapter 5: Writing"
---

Once you have generated the boilerplate template files for your module, you can start updating these to make your module work.
The boilerplate template files have lots of TODO comments to help give pointers how to write the module, as well example contents.

In this chapter we will go through each one, explaining what each section of each file is doing and why some bits are the way they are due to nf-core specifications.

## The `environment.yml` file

The first file that is generated is a conda `environment.yml` file.

If the `nf-core modules create` command was successfully able to find the tool within bioconda, it should already have specified the channel, tool name, and the specific version of the tool, as you can see in the example below.

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

The tool and version specified in this file will in most cases have a corresponding biocontainer (specified in `main.nf`).

You only need to modify this if your tool does not use tools from bioconda (e.g. a more generic tool that is only on `conda-forge`), and there is not an existing biocontainer image for the tool.
However, the nf-core guidelines strongly recommends conda support (as the most accessible software management system), so if your tool is not on bioconda or conda-forge, we strongly recommend you add this to the appropriate repository.

:::info{title="Behind the scenes" collapse}
nf-core uses a separate conda file rather than defining within the `main.nf` script file's conda directive, for [TODO] reason.
:::

## The `main.nf` file

This is the main file of the nf-core wrapper, which is where the Nextflow module code itself is defined.

When you first generate the file, some parts of the module will already be filled in for you, and there will be many TODO comments to help guide you through writing the module.

There are 9 main Nextflow process blocks of an nf-core module

- tag
- label
- conda
- container
- input
- output
- when
- script
- stub

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

The `label` block is automatically propogated by the option you specified during the boilerplate file generation.

```nextflow
    label 'process_single'
```

nf-core defines a preset set of labels with associated default resources, as specified [here](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L22-L54).
You should select from one of these standard labels that approximately matches the typical requirements of the tool you are creating the module for.

You should not use custom label names, as the default resource specifications for that label can be modified by a pipeline developer at a pipeline level (in nf-core pipelines, this happens in the `conf/base.config` file).
Furthermore a user can also override these within their own Nextflow custom config `process` block by specifying the module name

This flexibility at the pipeline level it is unnecessary to modify this away from the standard labels specified here, and provides the benefit of 'familiarity' between all modules and nf-core pipelines on how each default resource is defined.

### The `conda` block

The `conda` block tells Nextflow to use the associated `environment.yml` file to generate the required conda environment for this module, when a user has specified to use conda for running a pipeline.

```nextflow
    conda "${moduleDir}/environment.yml"
```

You should not modify this, except in the very rare exceptional cases where a given tool does not support conda (typically proprietary tools), in which case you may remove the line.

### The `container` block

The `container` block describes where to pull the tool's docker and singularity containers.

```nextflow
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0':
        'biocontainers/drep:3.5.0--pyhdfd78af_0' }"
```

In most cases the nf-core tool's module create command should have auto-detected based on the tool name the latest version and build of the given container, as with the conda `environment.yml` file.

You should not modify this, except in the very rare exceptional cases where a given tool does not have a bioconda recipe nor biocontainer (typically proprietary tools), in which case you may manually modify these lines to point to the corresponding container.

### The `input` block

The `input` block is where you will define what inputs the module will receive.

```nextflow
    input:
    tuple val(meta), path(bam)
```

The boilerplate template comes with an example of a common bioinformatic file format (`.bam`) and an associated `meta` map (the meta map may be missing if you did not specify this during the boilerplate template file creation).

You should edit this section of the `main.nf` to define the files your module will require.

The input block of an nf-core module should include all possible files via a `path()` declaration a given tool(/subcommand) can accept - regardless of whether these are mandatory or optional.
You should make sure to check the documentation of your tool to ensure you capture all of these!

In addition to files, mandatory tool argument/parameters (non-file) should also be specified here as a `val` input channel (in most cases without a meta map).
These are often required arguments without a default such as a `--mode` parameter.
Optional parameters can be specified via different mechanism described later, so only mandatory non-file arguments should be specified in input channels.

:::info{title="Behind the scenes" collapse}

We require _all_ files to be specified here because Nextflow processes _must_ receive files via this route to correctly 'stage' the files in a jobs working directory.

Furthermore, you cannot assume that your use case of the module will cover what other future pipeline developers may need to use it for (if you feel you only need a few files), thus you should try to be as comprehensive as possible.

While there is no 'native' Nextflow way to specify optional inputs, a pipeline developer can instead use an empty list `[]` to the module to indicate that file is not being provided.

We also require all mandatory non-file arguments to ensure that a particular module has all the minimum required information for a tool to run, so it can run 'out of the box' without any extra configuration (see `ext.args` below).
:::

The names of the inputs should generally follow the convention of using the file format as the variable name or the tool's argument name that is used in the command to refer to the input object.

If the tool has multiple input files, you should generally add these in their own input channels.
Each new channel can have their own `meta`, with each subsequent one being called `meta2`, `meta3`.
Generally we recommend all input channels should have a meta, except if you can guarantee that you would never have more than a single version of given input file (e.g., a configuration file that applies to every single sample/input data file).

In a few exceptional cases, multiple very closely related files can be specified 'inline' in a single channel, such as index files with the main file (e.g. `.bai` file with a `.bam` file), where the index file cannot be used by itself without the main file.

:::info{title="Behind the scenes" collapse}
The purposes of the 'one channel = one file' specification is to ensure as much consistency as possible for pipeline developers.

If we allowed high heterogeneity in number and ordering of files in a single input channel, it makes it harder for pipeline developers to know how to prepare input/output channels.

For a pipeline developer, they can use Nextflow's [`.multiMap()`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) operator to split their own multi-element channel into a separate 'sub-channels' to ensure the objects in the elements stay associated with each other even though being separated into channels.
:::

:::tip{title="Examples" collapse}

An example of an `input` block of a tool that requires a single input file and a mandatory argument:

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

The `output` block is where you will define all possible output files a tool can produce.

```nextflow
    output:
    path "versions.yml"           , emit: versions
```

The boilerplate template code by defaults includes a mandatory `versions.yml` file emission.
This is because all nf-core pipelines must report the versions of every tool used within the pipeline.

Above this line you can include the output definitions of all files your tool(/subcommand) will produce.
You should already define as many output files as possible (whether produced by default, or optionally), to make as much functionality available to future users as possible.
If your module uses a meta, you should emit all files with a meta.

Each output file type should have it's own output `emit` entry, e.g. one per file format.
In some cases if two mutually exclusive alternative file formats can be created by the tool, but serve the same purpose (such as two different indexing formats, for example `.bai` and `.csi` for `.bam files`), these can be captured in the same pattern and output channel.

:::info{title="Behind the scenes" collapse}
We often will recommend including as many input/output files as possible in the first creation of a tool.
This is because want to both make sure your nf-core module can be used by many people as possible, but also to minimise the number of channel structure changes/updates a module has to go through.

By being as comprehensive as possible during initial creation, we maintain stability and familiarity of how the module is used by pipeline developers and thus reduces the risk of breaking pipelines when they update modules.
:::

Generally, we recommend that where possible all output files should be compressed (if a tool can accept compressed inputs), thus you should ensure the path patterns take this into account.

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

You should not modify or remove this block in any case.

### The `script` block

The `script` block is the main part of the module, where the actual command that the module will perform within pipelines is defined.
You will do most modification of the `main.nf` here.

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

The boilerplate template already comes with an example samtools command for you in both the command error and within the `versions.yml` HEREDOC, which you will need to replace.

At the top of the block you have two standard variables that should not be removed.

- The`args` variable is how a pipeline developer can 'inject' optional parameters into the command itself.
  - The `ext.args` comes from a process scope of a nextflow configuration file - in nf-core pipelines this is in the `modules.config`
  - More information on this is described in the 'using in pipelines' chapter
- The `prefix` variable corresponds to the default output file basename that all output files should use
  - The default corresponds to the default `id` value of the primary input channel's metamap that all nf-core pipelines use.

Within the script section, you can replace the example with your own command, split across multiple lines with escaped backslashes for readability, specifying all the variables you need based on your input channels.

All commands must have the `$args` variable used, and if supported, the number of CPUs the command can use.

You can add additional variables after the `args` and `prefix` if you need to dynamically inject input files in different ways (e.g., in genomics if you have to use different flags for specifying single- or paired-end sequencing data files)

You must keep the `versions.yml` HEREDOC, but replace the command to generate a clean 'version' (i.e. just `1.2.3`, without leading `v` etc) of the version from the tools `--version` output.
We use standard UNIX command line tools for such a clean up, such as `cut` or `sed` for this purpose.
If your tool does not have way of emitting this information, you can use specify a dedicated variable that is interpreted within the HEREDOC.

While there is a relatively strict 'one module == one command/tool' rule, in some exceptional cases you may need to chain multiple tools (e.g. piping) or have a second command (e.g. compression).
In these cases you must define a second `args` variable with `args2`, `args3` and so on for each part of the pipe or tool.

It is important that you do not use any custom `meta` elements in modules!

:::warning

We do not provide any explicit examples here as there is a lot of variation across different tools.
We recommend you to just check multiple modules within the [nf-core/modules GitHub repository](https://github.com/nf-core/modules) to get further guidance.
:::

:::info{title="Behind the scenes" collapse}
Custom meta fields will not be standardised across tools, and add extra burden on pipeline developers who may want to encode the information in their metamaps by different names.

All optional/additional information for dynamic construction of a command away from the bare minimum execution should be left to the pipeline developer to inject via the `ext.args`
:::

### The `stub` block

The `stub` block is meant to just simulate the outputs emissions of a module during a 'dry run' of the pipeline.

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

This block should essentially copy the `script` block, but instead executing the actual command, you just create empty files with the same names as actual output files the tool would create using one or more `touch` commands.
If your tool outputs gzipped files, you should pipe the output of an empty echo command into `gzip`,

```bash
echo "" | gzip > ${prefix}.txt.gz
```

You should create files that can be picked up by every `output` channel.

If you get [Nextflow language server](https://www.nextflow.io/docs/latest/vscode.html) errors, you can either copy the command from the `script` block and wrap it in an `echo`, or you may simplify remove the `args` definition.

It is important to still copy/update the `versions.yml` HEREDOC command, as this will still be executed during the dry run.

## The `meta.yaml` file

Another nf-core specific component of a module is the `meta.yaml` file.

The purpose of this file is to better document all components of the module, with additional descriptions, keywords, and links to the tool's original resources.
Much of this information is used to help make the modules searchable on the nf-core website [modules page](https://nf-co.re/modules), make the wrapper findable on external websites, and to prepare to move towards automate linkage between modules.

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

This section provides metadata about the tool and module itself to make it more findable e.g. on the nf-core website.

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

Some parts of this section are automatically filled in for you based on the name of the module, and information already recorded in metadata of the bioconda recipe that is used in the previously described `environment.yaml`

In most cases, you only have to update the keywords, and the module description.

For the module description (second line after the `name`), you need to describe the purpose of the _module_ in a short sentence.

For the keywords, you need to have a minimum of 3 - all lower case - and should include things such as the tool's name, subcommand, words to describe the 'action' the module is performing such as 'sort' or 'filter', but can also be words describing the research field the module used - such as genomics, metagenomics, earth sciences, and the file formats the action is used on.
Ultimately, the more keywords the better.

Sometimes you may need to fill in the URLs to the tool's own source code and documentation.
You can just copy and paste these in after a quick Google search.
It is OK to leave these empty or duplicate URLs if the tool doesn't have e.g. a dedicated documentation page.

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

In the default example added by the boilerplate template, you see a single input channel described - one with a metamap and a single path element.

Here you need to update the the input channel name, the description, pattern, and ontologies.

The type follows a fixed set of [categories](https://nf-co.re/docs/guidelines/components/modules#input-and-output-channel-types), such as `file`, `integer`, `boolean`.
The description should be descriptive in what they contain or how they should be prepared, e.g. 'A TSV file containing 5 columns generated by XYZ', not just a `TSV`.
The pattern should match the corresponding input expected by the tool itself.
The ontologies section should provide links to respective entry matching the input type within the the controlled vocabulary of the (typically) [EDAM ontology](https://www.ebi.ac.uk/ols4/ontologies/edam).
It is up to you how specific to be - whether extremely specific (fastq), or generic (text file).

:::tip
You can use the `nf-core modules lint` command (described later) to automatically pick up all the input channels and put placeholders elements in for you
:::

Note that you will have to make an entry for each of your input channels
If you copy and paste, make sure to edit each section to be unique for each of your input channels.
When each of the additional meta files have their own `meta` maps, the name should be updated to `meta2` etc accordingly.

### Output channel descriptions

This section follows the same concept as the input section.

```
## TODO nf-core: Add a description of all of the variables used as output
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
You can use the `nf-core modules lint` command (described later) to automatically pick up all the input channels and put placeholders elements in for you
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

## Summary

In this chapter we have gone through each file in fine detail, describing how each part of each file should be filled out.
This tutorial cannot cover every single point of the specifications, nor cover every edge case.
We therefore highly recommend to always refer to the [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules) when writing a module.

In the next chapter we will describe the other half of an nf-core modules - the files used for unit testing of the modules.
