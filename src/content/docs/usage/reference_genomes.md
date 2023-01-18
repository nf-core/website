---
title: Reference Genomes
subtitle: How reference genomes are handled in nf-core
---

Many nf-core pipelines need a reference genome for alignment, annotation, or similar purposes.

## Illumina iGenomes

To make the use of reference genomes easier, Illumina has developed a centralised resource called [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).
The most commonly used reference genome files are organised in a consistent structure for multiple genomes.

We have uploaded a copy of iGenomes onto AWS S3 and nf-core pipelines are configured to use this by default.
AWS-iGenomes is hosted by Amazon as part of the [Registry of Open Data](https://registry.opendata.aws/aws-igenomes/) and are free to use. For more information about the AWS iGenomes, see [https://ewels.github.io/AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/).

All AWS-iGenome paths are specified in pipelines that support them.
By default, the pipeline will automatically download the required reference files when you run the pipeline and supply an appropriate genome key (eg. `--genome GRCh38`).

Downloading reference genome files takes time and bandwidth so, if possible, we recommend making a local copy of the iGenomes resource.

To use a **local** version of iGenomes, the variable `params.igenomes_base` must be set to the path of the local igenome folder. To do so, you can use either:

- :white_check_mark: a configuration profile;
- :white_check_mark: the `--igenomes_base` parameter in your execution command;
- :x: setting the `params.igenomes_base` path in a custom config file will **not** work. Parameters should be provided within a `params file` with the option `-params-file` instead.

Additional information on how to set a local igenomes base can be found [here](troubleshooting.md#using-a-local-version-of-igenomes).

## Adding paths to a config file

If you are unable to use the AWS-iGenomes references, you can still supply reference genome paths on the command line via the pipeline's parameters e.g. `--fasta` or `--gtf`.

If you are using the same references repeatedly, it can be more convenient to save these paths in a nextflow config file.
Pipelines that support AWS-iGenomes can also be configured to support custom genome IDs and paths.

To use this system, add paths to your config file using the following template:

```nextflow
params {
  genomes {
    'YOUR-ID' {
      fasta  = '/path/to/data/genome.fa'
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```

You can add as many genomes as you like as long as they have unique IDs.
References are used with the command line option `--genome YOUR-ID`.

Read the [Nextflow configuration documentation](configuration.md) for more information about custom config files.

## Using Refgenie for genome management

You can also use the reference genome manager [Refgenie](http://refgenie.databio.org/en/latest/overview/) with nf-core pipelines.

1. Install and initialize refgenie following the official [documentation](http://refgenie.databio.org/en/latest/install/).

A file required by nf-core containing refgenie genome assets will be automatically created at `~/.nextflow/nf-core/refgenie_genomes.config`. An `includeConfig` statement to this new config file will be appended to the `~/.nextflow/config` file. This file is automatically loaded by Nextflow during every pipeline run (see the [Nextflow documentation](https://nextflow.io/docs/latest/config.html)).

To use a new reference genome or asset, fetch it via normal refgenie usage (`refgenie pull`) - the nf-core plugin will automatically update the `refgenie_genomes.config` configuration file.
This file should never be edited manually, as it is overwritten during each refgenie command.

2. Pull all the reference assets that you may need to run the pipeline.

```bash
refgenie pull t7/fasta
refgenie pull t7/bowtie2_index
```

Asset paths are automatically added to `~/.nextflow/nf-core/refgenie_genomes.config`, included in `~/.nextflow/config` and available to every pipeline run.

The file format mimics the `igenomes.config` file that comes with many nf-core pipelines:

```nextflow
// This is a read-only config file managed by refgenie. Manual changes to this file will be overwritten
// To make changes here, use refgenie to update the reference genome data
params {
  genomes {
    't7' {
      bowtie2_index        = "<path to refgenie genomes>/alias/t7/bowtie2_index/default/t7"
      fasta                = "<path to refgenie genomes>/alias/t7/fasta/default/t7.fa"
    }
  }
}
```

Here, the genome _key_ that you'll use to launch the pipeline is `t7`.

> NOTE: You can also use [custom assets](http://refgenie.databio.org/en/latest/custom_assets/).

3. Run your pipeline, specifying the required genome.

```bash
nextflow run nf-core/<PIPELINENAME> --genome t7  # ..rest of normal pipeline flags..
```

Please refer to [Refgenie documentation](http://refgenie.databio.org/en/latest/) for further information.
