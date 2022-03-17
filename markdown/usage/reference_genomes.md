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

To use a local version of iGenomes, the variable `params.igenomes_base` must be set to the path of the local igenome folder using either a configuration profile or by using `--igenomes_base` in your execution command.
Setting the `params.igenomes_base` path in a custom config file will **not** work.
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
