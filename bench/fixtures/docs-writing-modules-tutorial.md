---
title: "Chapter 5: Writing your module"
subtitle: "How to fill in the module files"
shortTitle: "Chapter 5: Writing modules"
---

Once you have generated the boilerplate template files, you can update them to make your module function.
The boilerplate contains many `TODO` comments and example content to guide you.

This chapter walks through each file, explaining what each section does and why, as defined by the nf-core specifications.

## The `environment.yml` file

The `environment.yml` file is a Conda environment specification.

If `nf-core modules create` found your tool on Bioconda, the file already contains the channel, tool name, and version:

```yaml
---
# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
channels:
  - conda-forge
  - bioconda
dependencies:
  - "bioconda::drep=3.5.0"
```

You rarely need to modify this file. The tool and version specified here will usually have a matching biocontainer in the `main.nf` file.

Only modify `environment.yml` when your tool is not on Bioconda and no recipe exists for it.
nf-core guidelines strongly recommend Conda support, as Conda is the most accessible software management system.
If your tool is not on Bioconda or [conda-forge](https://conda-forge.org/), add it to the appropriate repository.
Bioconda provides [tutorials](https://bioconda.github.io/tutorials/2024-adding-bioinformatic-software-to-bioconda.html) on their website.

:::info{title="Behind the scenes" collapse}
nf-core uses a separate Conda file rather than the `main.nf` conda directive for two reasons:

- It enables [automated container building with Seqera's wave infrastructure](https://nf-co.re/blog/2024/seqera-containers-part-1).
- It makes multi-tool environments easier to read and manage than multiple conda declarations on a single line.
  :::

## The `main.nf` file

The `main.nf` file contains the Nextflow module code.
When generated, parts of the module are pre-filled and `TODO` comments mark where you need to make changes.

An nf-core module has 9 main Nextflow process blocks:

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

The boilerplate `TODO` comments have been removed for readability.

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

Each block is described below.

### The `tag` block

The `tag` block propagates standard metadata from each input channel's [meta map](https://nf-co.re/docs/developing/components/meta-map).

```nextflow
tag "$meta.id"
```

You generally do not need to modify this.
All nf-core pipelines assume an `id` element in each meta map, and pipeline developers can customise it.

### The `label` block

The `label` block is pre-filled with the label you selected during boilerplate generation.

```nextflow
label 'process_single'
```

Key points:

- nf-core defines a preset set of [standard labels with default resources](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L22-L54).
- Select the label that best matches your tool's typical resource needs.
- Do not use custom labels. Pipeline developers can tune default resources in the pipeline's `conf/base.config`, and users can override per-module resources using `withName:` in a custom config.
- Standard labels provide modules and pipelines a consistent set of default resources.

### The `conda` block

The `conda` block tells Nextflow to use the `environment.yml` file when a pipeline runs with Conda.

```nextflow
conda "${moduleDir}/environment.yml"
```

Do not modify this.
You may only remove this line when a tool that does not support Conda (typically proprietary tools).

### The `container` block

The `container` block tells Nextflow where to pull the tool's Docker and Singularity containers from.

```nextflow
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0':
    'biocontainers/drep:3.5.0--pyhdfd78af_0' }"
```

The `nf-core modules create` command auto-detects the container name, version, and build from Bioconda and biocontainers.

Do not modify this unless the tool has no Bioconda recipe or biocontainer (typically proprietary tools).
In that case, update the links to the relevant container.

### The `input` block

The `input` block defines what the module receives.

```nextflow
input:
tuple val(meta), path(bam)
```

The boilerplate includes an example BAM input and a [meta map](https://nf-co.re/docs/developing/components/meta-map).
Edit this section to declare every file your module needs.

Guidelines:

- Declare **all possible files** with a `path()` entry, both mandatory and optional. Check the tool's documentation to capture them all.
- Declare mandatory non-file arguments as `val` input channels, usually without a meta map.
- Name channels after the file format or the tool's command-line argument.
- Use one channel per input file.
- Each extra input channel gets its own meta map, named `meta2`, `meta3`, and so on.
- Every input channel should have a meta, unless you can guarantee only one version of that file exists across all samples (for example, a shared configuration file).
- Closely related files can share a tuple channel when neither can be used alone, such as a `.bai` index with a `.bam` file.

```nextflow
input:
tuple val(meta), path(bam)
val mode
```

```nextflow
input:
tuple val(meta), path(bam), path(bai)
```

:::info{title="Behind the scenes" collapse}
Nextflow processes must receive files through input channels to stage them correctly in the job's working directory.
Declaring every possible input keeps the module reusable across pipelines, and minimises breakages for other developers.

Nextflow has no native optional input syntax. Pipeline developers pass an empty list `[]` to indicate an input is not provided.

Mandatory non-file arguments belong in the input block so the module can run without extra configuration.
Optional parameters are handled through `ext.args`, described below.

The one-channel-one-file rule keeps channel structures consistent.
See the nf-core [module specifications](../../../specifications/components/overview) for exceptions.
Pipeline developers can use Nextflow's [`.multiMap()`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) operator to split multi-element channels into sub-channels while keeping elements associated.
:::

:::tip{title="Examples" collapse}

Single input file with a mandatory argument:

```nextflow
input:
tuple val(meta), path(fasta)
val(output_format)
```

Multiple input files:

```nextflow
input:
tuple val(meta), path(reads)
tuple val(meta2), path(fasta)
tuple val(meta3), path(fasta_index)
```

Input with a closely associated index file:

```nextflow
input:
tuple val(meta), path(input), path(intervals)
path  fasta
```

:::

### The `output` block

The `output:` block declares every file the tool can produce.

```nextflow
output:
tuple val("${task.process}"), val('<tool1>'), eval('tool1 --version'), emit: versions_tool1, topic: versions
```

The boilerplate always includes a mandatory 'topics' version emission channel because nf-core pipelines must report every tool version.

Guidelines:

- Declare as many outputs as possible during first creation, whether default or optional, to maximise reuse and minimise future channel changes.
- If the module uses a meta map, emit all files with it.
- Specify a command to get a 'clean' version number (e.g. just `1.2`, not `samtools v1.2 (2026-01-02)`) in the `eval()` section of the 'versions' topic channel.
- Give each file type its own `emit` entry, usually one per format.
- Combine two mutually exclusive formats that serve the same purpose (such as `.bai` and `.csi` indexes for `.bam`) in one channel.
- Compress output files where the tool supports compressed input, and reflect this in the path pattern.
- Name `emit:` channels after the file format or suffix.
- Mark variable outputs with Nextflow's `optional: true`.

```nextflow
output:
tuple val(meta), path("${prefix}.bam"),  emit: bam,  optional: true
tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
tuple val(meta), path("${prefix}.sam"),  emit: sam,  optional: true
tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | cut -d ' ' -f 2'), emit: versions_samtools, topic: versions
```

```nextflow
output:
tuple val(meta), path("*.{vcf,vcf.gz}"),  emit: vcf
```

:::info{title="Behind the scenes" collapse}
Comprehensive inputs and outputs during initial creation keep modules reusable and reduce breaking changes. Pipeline developers rely on a stable channel structure.
:::

:::tip{title="Examples" collapse}
A module producing multiple files, some optional:

```bash
output:
tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
tuple val(meta), path('*.json')           , emit: json
tuple val(meta), path('*.html')           , emit: html
tuple val(meta), path('*.log')            , emit: log
tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | sed -e "s/fastp //g"'), emit: versions_fastp, topic: versions
```

A module producing three mutually exclusive primary outputs:

```nextflow
output:
tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
tuple val(meta), path("*.csi")                    , emit: csi, optional: true
tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools
```

A module producing multiple files that all support compressed variants:

```nextflow
output:
tuple val(meta), path('*.{blast,blast.gz}'), optional: true, emit: blast
tuple val(meta), path('*.{xml,xml.gz}')    , optional: true, emit: xml
tuple val(meta), path('*.{txt,txt.gz}')    , optional: true, emit: txt
tuple val(meta), path('*.{daa,daa.gz}')    , optional: true, emit: daa
tuple val(meta), path('*.{sam,sam.gz}')    , optional: true, emit: sam
tuple val(meta), path('*.{tsv,tsv.gz}')    , optional: true, emit: tsv
tuple val(meta), path('*.{paf,paf.gz}')    , optional: true, emit: paf
tuple val("${task.process}"), val('diamond'), eval('diamond --version 2>&1 | tail -n 1 | sed "s/^diamond version //"'), emit: versions_diamond, topic: versions
```

:::

### The `when` block

The `when` block lets pipeline developers dynamically activate or deactivate the module.

```nextflow
when:
task.ext.when == null || task.ext.when
```

Do not modify or remove this block.

### The `script` block

The `script` block defines the command the module runs.
You will edit this section most.

```nextflow
script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
samtools \\
    sort \\
    ${args} \\
    -@ $task.cpus \\
    -o ${prefix}.bam \\
    -T $prefix \\
    $bam
"""
```

The boilerplate includes an example samtools command.
Replace with the command for your tool.

Two standard variables appear at the top of the block. Do not remove them:

- `args` — how pipeline developers inject optional parameters into the command. The value comes from `ext.args` in the process scope of a Nextflow configuration file (defined in `modules.config` for nf-core pipelines). See the ["Using in pipelines"](./8-using) chapter for details.
- `prefix` — the default output file basename. It defaults to the `id` value of the primary input channel's meta map.

Replace the example command with your own command, split across multiple lines with escaped backslashes for readability.
Reference all required input variables.

Every command must use the `$args` variable, and the CPU count where supported.

You can add extra variables after `args` and `prefix` to inject input files dynamically. For example, to handle single- or paired-end sequencing data:

```nextflow
script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
def single_end  = meta.single_end ? "--single" : ""
"""
elprep merge \\
       input/ \\
       output/${prefix}.bam \\
       ${args} \\
       ${single_end} \\
...
"""
```

If you must chain tools (for example piping or compression), define additional `args2`, `args3`, and so on, one per command in the pipe.

Do not use custom `meta` elements in modules.

:::info{title="Behind the scenes" collapse}
Custom meta fields are not standardised across tools and create extra work for pipeline developers.
All optional parameters for dynamic command construction should come from `ext.args`.
:::

:::warning
This chapter does not include explicit `script` block examples because tool commands vary widely. Browse the [nf-core/modules GitHub repository](https://github.com/nf-core/modules) for reference modules.
:::

### The `stub` block

The `stub` block simulates the module's output during a `-dry-run` of a pipeline.

```nextflow
  stub:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  touch ${prefix}.bam
  """
```

The stub block should mirror the `script` block, except instead of running the tool, it uses `touch` to create empty files matching every output name.
For gzipped outputs, pipe an empty echo into `gzip`:

```nextflow
stub:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
echo "" | gzip > ${prefix}.txt.gz
"""
```

Create files that every `output` channel can pick up.

If you edit in VS Code and see [Nextflow language server](https://www.nextflow.io/docs/latest/vscode.html) errors, either copy the command from the `script` block and wrap it in `echo`, or remove the `args` definition:

```nextflow
stub:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"
"""
echo "$args"
echo "" | gzip > ${prefix}.txt.gz
"""
```

## The `meta.yml` file

The `meta.yml` file documents the module with descriptions, keywords, and links to the tool's resources.
This information powers the searchable [modules page](https://nf-co.re/modules), improves discoverability on external databases, and supports future automated linkage between modules.

The file has four main sections:

- Tool and module description.
- Input channel descriptions.
- Output channel descriptions.
- Contributors.

:::info{title="Click here to see full 'raw' file example" collapse}

The boilerplate `TODO` comments have been removed for readability.

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

  versions_drep:
    - - ${task.process}:
          type: string
          description: The name of the process
      - drep:
          type: string
          description: The name of the tool
      - drep version 2>&1 | sed 's/^.*drep v//':
          type: eval
          description: The expression to obtain the version of the tool
topics:
  versions:
    - - ${task.process}:
          type: string
          description: The name of the process
      - drep:
          type: string
          description: The name of the tool
      - drep version 2>&1 | sed 's/^.*drep v//':
          type: eval
          description: The expression to obtain the version of the tool

authors:
  - "@jfy133"
maintainers:
  - "@jfy133"
```

:::

### Tool and module description

This section holds metadata for finding the module on the nf-core website.

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

Much of this section is pre-filled from the module name and the Bioconda recipe referenced in `environment.yml`.

What you need to update:

- **`description`** — a short sentence describing the module's purpose.
- **`keywords`** — a minimum of three, all lowercase. Include the tool name, subcommand, the action the module performs (such as `sort` or `filter`), the research field (such as `genomics` or `metagenomics`), and the file formats involved. More keywords improve discoverability.
- **URLs** — fill in source code and documentation links. Leave empty or duplicate where the tool has no dedicated documentation page.
- **Multiple tools** — if your module uses more than one tool (for example the tool and `gzip`), add a description block for each.

### Input channel descriptions

This section describes every input channel.

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

The boilerplate template shows a single input channel with a meta map and a path element.

For each channel, update:

- **Name** — match the input channel variable.
- **`type`** — one of the [fixed categories](https://nf-co.re/docs/specifications/components/modules/input-output-options) (`file`, `integer`, `boolean`, and so on).
- **`description`** — describe the contents or preparation, for example "A TSV file containing 5 columns generated by XYZ", not just "TSV".
- **`pattern`** — match the input the tool expects.
- **`ontologies`** — link to matching entries in a controlled vocabulary, typically the [EDAM ontology](https://www.ebi.ac.uk/ols4/ontologies/edam). Choose the level of specificity that suits your input.

Add one entry per input channel. When a channel has its own meta map, rename it to `meta2`, `meta3`, and so on.

:::tip
Run `nf-core modules lint` to auto-populate placeholders for every input and output channel.
:::

:::tip{title="Examples" collapse}
A complete `meta.yml` with ontologies and multiple meta entries (see below for `outputs:`):

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
  versions_bwamem2:
    - - ${task.process}:
          type: string
          description: The name of the process
      - bwamem2:
          type: string
          description: The name of the tool
      - bwa-mem2 version | grep -o -E "[0-9]+(\.[0-9]+)+":
          type: eval
          description: The expression to obtain the version of the tool
  versions_samtools:
    - - ${task.process}:
          type: string
          description: The name of the process
      - samtools:
          type: string
          description: The name of the tool
      - samtools version | sed '1!d;s/.* //':
          type: eval
          description: The expression to obtain the version of the tool
topics:
  versions:
    - - ${task.process}:
          type: string
          description: The name of the process
      - bwamem2:
          type: string
          description: The name of the tool
      - bwa-mem2 version | grep -o -E "[0-9]+(\.[0-9]+)+":
          type: eval
          description: The expression to obtain the version of the tool
    - - ${task.process}:
          type: string
          description: The name of the process
      - samtools:
          type: string
          description: The name of the tool
      - samtools version | sed '1!d;s/.* //':
          type: eval
          description: The expression to obtain the version of the tool
authors:
  - "@maxulysse"
  - "@matthdsm"
maintainers:
  - "@maxulysse"
  - "@matthdsm"
```

:::

### Output channel descriptions

This section follows the same pattern as inputs.

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

The only difference: the output file key (not the channel name) must match the pattern used in the output channel in `main.nf`.

Update the `type`, `description`, `pattern`, and `ontologies` as for inputs.

### Contributors

This section lists the original author and any subsequent contributors.

```yaml
authors:
  - "@jfy133"
maintainers:
  - "@jfy133"
```

The boilerplate populates this automatically. You do not need to modify it.

If you later update a module originally written by someone else, add yourself to `maintainers:`.

## Further reading

This chapter cannot cover every specification detail or edge case. Always refer to the [nf-core module specifications](https://nf-co.re/docs/specifications/components/modules/general) when writing a module.

The [next chapter](./6-testing) covers unit testing the module with nf-test.
