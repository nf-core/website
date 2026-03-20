---
title: Meta maps and sample metadata
subtitle: Using sample metadata stored in meta maps
shortTitle: Meta maps
---

nf-core pipelines use 'meta maps' to carry sample metadata alongside files throughout a pipeline run.

They are commonly used for recording IDs or names, tracking pipeline-generated metadata, grouping or splitting operations, and to specify per-sample tool arguments.
For a detailed practical tutorial on how to create and use meta maps within pipeline code, refer to the dedicated Nextflow training tutorial [side quest](https://training.nextflow.io/latest/side_quests/metadata/#12-pick-out-specific-fields-with-map).

## Structure

A meta map in nf-core pipelines are groovy a [unordered key-value list](https://www.tutorialspoint.com/groovy/groovy_maps.htm).

```nextflow
[id: 'test', single_end: false]
```

A meta map sits a tuple within a Nextflow channel object, next to the one or more files the metadata is describing.

```nextflow
[ [id: 'test', single_end: false], sample1_R1.fastq.gz ]
```

```nextflow
[ [id: 'test', single_end: true], [ sample1_R1.fastq.gz, sample1_R2.fastq.gz], reference.fasta ]
```

## Key names

nf-core developers may define and within a pipelines or local subworkflows any name for a meta map key, and record any metadata they require for the execution of the pipeline.

nf-core only defines 2 'standard' meta map keys.

| key               | purpose                                                                                          |
| ----------------- | ------------------------------------------------------------------------------------------------ |
| `meta.id`         | recording unique file identifiers associated with a file (e.g. 'sample' names in bioinformatics) |
| `meta.single_end` | genomic sequencing pipelines handling paired-end sequencing data.                                |

There are no other standard or required key names that nf-core developers need to use.

:::note{collapse title="Note on meta.single_end"}
The standard `meta.single_end` key is a relic from the early development of meta maps.
The key became common across many pipelines due to the presence of the key as an example in the nf-core template.
No other 'standard' meta keys will be officially defined in the future.
This is to provide maximum flexibility to pipeline developers.
:::

## Usage in nf-core components

### Modules

The two standard meta map keys (`id` and `single_end`) are the [only keys allowed](../../specifications/components/modules/general#types-of-meta-fields) to be explicitly referred to in an nf-core/module.

nf-core/modules refer to meta maps in process `input:` blocks via an entry within a tuple.

```groovy
input:
tuple val(meta), path(reads)
```

These keys are used to define a unique process identifier with the Nextflow `tag` directive:

```nextflow title="main.nf"
process FASTQC {
    tag "${meta.id}"
...
```

And define the default file name `prefix`:

```nextflow title="main.nf"
...
    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
...
```

All other usage of meta map keys within a module must come via the `ext.args` variable, and specified in `modules.config`.

### Subworkflows

No meta data keys are to be assumed to be present in input channels other than the standard [key names], similarly to [modules](#modules).

In contrast to modules, nf-core/subworkflows are allowed to generate new meta map keys, that can be optionally emitted at the end of the subworkflow.

New meta map keys defined during the subworkflow must be documented in the associated `meta.yml` of the nf-core/subworkflow.

```yaml title="meta.yaml"

...
output:
  - bam:
      description: |
        Structure: [ val(meta), path(bam) ]
        Sorted BAM/CRAM/SAM file
        Note: the subworkflow adds a new meta ID `meta.id_index` that _must_
        be used in `prefix` to ensure unique file names. See the associated
        nextflow.config file.
...
```

## Usage in pipelines

Pipeline developers may define and use meta maps throughout their code as necessary.

Official nf-core/modules or nf-core/subworkflows installed within the pipelines should not be modified to add custom meta map keys.

Usage of meta maps information within modules are restricted to `modules.conf`.
The information encoded by the keys can then be passed by `ext.args` and other process directives.

```nextflow title='modules.config'
process {
    withName: MALT_RUN {
        tag        = { "${meta.db_name}|${meta.id}" }
        ext.args   = { "${meta.db_params} -m ${params.malt_mode}" }
        ext.prefix = { "${meta.db_name}" }
        publishDir = [
            path: { "${params.outdir}/malt/${meta.db_name}/" },
            mode: params.publish_dir_mode,
            pattern: '*.{rma6,log,sam}',
        ]
    }
}
```

## Useful links

- **[Nextflow training tutorial](https://training.nextflow.io/latest/side_quests/metadata/#12-pick-out-specific-fields-with-map)**: A step-by-step guided tutorial for writing and using meta maps
- **[Rob Syme's nf-core/bytesize talk](https://www.youtube.com/watch?v=A357C-ux6Dw&t=597s)**: A youtube video on how to safely modify Groovy meta map objects within Nextflow
