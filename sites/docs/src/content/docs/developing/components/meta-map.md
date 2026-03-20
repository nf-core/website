---
title: Meta maps and sample metadata
subtitle: Carrying sample metadata through pipelines with meta maps
shortTitle: Meta maps
---

<!-- TODO Initially started with a concept, but this is getting bigger, can't decide whether to switch to remove detail or switch to reference  -->

nf-core pipelines use 'meta maps' to carry sample metadata alongside files throughout a pipeline run.
They are commonly used for recording IDs or names, tracking pipeline-generated metadata, grouping or splitting operations, and to specify per-sample tool arguments.

## Structure

A meta map in nf-core pipelines are groovy a [unordered key-value list](https://www.tutorialspoint.com/groovy/groovy_maps.htm).

```groovy
[id: 'test', single_end: false]
```

A meta map sits a tuple within a Nextflow channel object, next to the one or more files the metadata is describing.

```groovy
[ [id: 'test', single_end: false], sample1_R1.fastq.gz ]
```

```groovy
[ [id: 'test', single_end: true], [ sample1_R1.fastq.gz, sample1_R2.fastq.gz], reference.fasta ]
```

## Key names

nf-core developers may define and within a pipelines or local subworkflows any name for a meta map key, and record any metadata they require for the execution of the pipeline.

The only pre-defined key names are `meta.id` for recording unique file identifiers, and `meta.single_end` for genomic sequencing pipelines handling paired-end sequencing data.
There are no other standard or required key names that nf-core developers need to use.

:::note{collapse title="Note on meta.single_end"}
The standard `meta.single_end` key is a relic from the early development of meta maps.
The key became common across many pipelines due to the presence of the key as an example in the nf-core template.
No other 'standard' meta keys will be officially defined in the future.
This is to provide maximium flexibility to pipeline developers.
:::

## Metamap Creation

Basic definition of a meta map can be with a list.

```groovy
def reads = [ [id: 'test', single_end: false], file(sample1_R1.fastq.gz) ]
```

In most cases however, meta maps will be created and updated within Nextflow channels.

Generation of a meta map within a channel occurs through using the `.map{}` operator.

```groovy
ch_fastq = Channel.fromPath('sample1_R1.fastq.gz')
  .map {
    fastq ->
      [ [id: fastq.getName() ], fastq ]
  }
```

## Usage in nf-core pipelines

Meta map key names are defined at a pipeline level in nf-core pipelines.

<!-- TODO

- In modules.config
- Modifying (remember Rob's )
- Combining
- Grouping

-->

## Usage in nf-core modules

nf-core modules refer to meta maps in process `input:` blocks via a tuple.

```groovy
input:
tuple val(meta), path(reads)
```

The two standard meta map keys (`id` and `single_end`) are the [only keys allowed](../../specifications/components/modules/general.md#types-of-meta-fields) to be explicitly referred to in an nf-core/module.

The keys are used to define a unique process identifier with the Nextflow `tag` directive, and define the default file name `prefix`.

All other usage of meta map keys within a module must come via the `ext.args` variable, and specified in `modules.config`.

## XXX

<!-- TODO: Rewrite this page, the content felt out of date and hard to follow

This page: https://nf-co.re/docs/contributing/components/meta_map#generating-a-meta-map-from-file-pairs

- [x] What is a metamap
- [x] What it is used for
- ?Common patterns?
  - [x] Creating
  - [ ] Grouping
  - [ ] Combining
  - [ ] Modifying
- For a tutorial, see nextflow training

-->

## Useful links

- **[Nextflow training tutorial](https://training.nextflow.io/latest/side_quests/metadata/#12-pick-out-specific-fields-with-map)**: A step-by-step guided tutorial for writing and using meta maps
- **[Rob Syme's nf-core/bytesize talk](https://www.youtube.com/watch?v=A357C-ux6Dw&t=597s)**: A youtube video on how to safely modify Groovy meta map objects within Nextflow
