---
title: Reference genomes
subtitle: Learn about nf-core reference genomes
shortTitle: Reference genomes
---

Many nf-core pipelines use reference genomes for alignment, annotation, and similar tasks.
This page describes available approaches for managing reference genomes.

There are three main ways to use reference genomes with nf-core pipelines:

- [Local copies of genomes](#local-copies-of-genomes): user downloaded and self-managed
- [AWS iGenomes](#aws-igenomes): Illumina-hosted pre-build reference genomes and indices
- [Refgenie](#refgenie): programmatic genome asset management tool

## Local copies of genomes

Most genomics nf-core pipelines can start from just a FASTA and GTF file and create downstream reference assets (genome indices, interval files, etc.) as part of pipeline execution.

Using GRCh38 as an example:

1. Download the latest files:

   ```bash
   #!/bin/bash

   VERSION=108
   wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
   wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz
   ```

2. Run pipeline with `--save_reference` to generate indices:

   ```bash
   nextflow run \
       nf-core/rnaseq \
       --input samplesheet.csv \
       --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
       --gtf Homo_sapiens.GRCh38.108.gtf.gz \
       --save_reference
   ```

   :::note
   The pipeline will generate and save reference assets. For example, the STAR index will be stored in `<results_dir>/genome/index/star`.
   :::

3. Move generated assets to a central, persistent storage location for re-use in future runs.
4. Use pre-generated indices in future runs.

   ```bash
   nextflow run \
       nf-core/rnaseq \
       --input samplesheet.csv \
       --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
       --gtf Homo_sapiens.GRCh38.108.gtf.gz \
       --star_index </path/to/moved/star/directory/> \
       --gene_bed </path/to/moved/genes.bed>
   ```

## AWS iGenomes

AWS iGenomes is Illumina's centralized resource that organizes commonly used reference genome and pre-built index files in a consistent structure for multiple genomes.
It provides the following benefits:

- Hosted on AWS S3 through the [Registry of Open Data](https://registry.opendata.aws/aws-igenomes/)
- Free to access and download
- Maintained by nf-core as a copy of the original Illumina resource

See the [AWS iGenomes documentation](https://ewels.github.io/AWS-iGenomes) for more information.

:::warning{title="Outdated annotations"}
Transcriptome and GTF files in iGenomes are significantly outdated.
For example, human annotations are from Ensembl release 75, while current release is 108+.
Consider using custom genomes for current annotations.
:::

:::warning{title="GRCh38 assembly issues"}
GRCh38 in iGenomes comes from NCBI instead of Ensembl, not the masked Ensembl assembly.
This can cause pipeline issues in some cases. See [nf-core/rnaseq issue #460](https://github.com/nf-core/rnaseq/issues/460) for details.
For GRCh38 with masked Ensembl assembly, use [Custom genomes](#custom-genomes).
:::

### Use remote AWS iGenomes

To use remote AWS iGenomes in supported nf-core pipelines, supply the `--genome` flag to your pipeline (e.g., `--genome GRCh37`).
On execution the pipeline will then:

1. Automatically download required reference files.
2. Auto-populated reference genome parameters from `conf/igenomes.config`.
   - Parameters like FASTA, GTF, and index paths are set automatically.
3. Download only what it requires for that specific workflow.

:::tip
Downloading reference genome files takes time and bandwidth.
We recommend using a local copy when possible. See [Use local AWS iGenomes](#use-local-aws-igenomes) for more information.
:::

### Use local AWS iGenomes

To use local AWS iGenomes:

1. [Download](https://github.com/ewels/AWS-iGenomes?tab=readme-ov-file#download-script) the iGenomes reference files you need to a local directory.
2. Set `--igenomes_base` to your local iGenomes directory path.

   :::warning
   This directory structure must reflect the structure defined in [`conf/igenomes.config`](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/igenomes.config).
   :::

3. Pipeline will use local files instead of downloading from AWS.

### Check annotation versions

To check the version of annotations used by AWS iGenomes:

1. Download the README file from the iGenomes S3 bucket using AWS CLI:

   ```bash
   aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt
   ```

2. View the README to see annotation details:

   ```bash
   cat README.txt
   ```

   Example output:

   ```console
   The contents of the annotation directories were downloaded from Ensembl on: July 17, 2015.

   Gene annotation files were downloaded from Ensembl release 75. SmallRNA annotation files were downloaded from miRBase release 21.
   ```

   This confirms the annotations are from Ensembl release 75 (July 2015), which is significantly outdated.

## Refgenie

Refgenie provides programmatic genome asset management as an alternative to manual file handling.

Refgenie allows you to:

- Automate asset tracking
- Manage genomes programmatically
- Use a consistent configuration format
- Update and version assets easily

### Use Refgenie

To use Refgenie:

1. Install Refgenie following the [official documentation](http://refgenie.databio.org/).
2. Initialize Refgenie.

   :::note
   Refgenie creates `~/.nextflow/nf-core/refgenie_genomes.config` and appends an `includeConfig` statement to `~/.nextflow/config` that references this file.
   :::

3. Pull required genome assets. For example:

   ```bash
   refgenie pull t7/fasta
   refgenie pull t7/bowtie2_index
   ```

   Asset paths are automatically added to `~/.nextflow/nf-core/refgenie_genomes.config`.
For example:

   ```groovy title="refgenie_genomes.config"
   // This is a read-only config file managed by refgenie. Manual changes to this file will be overwritten.
   // To make changes here, use refgenie to update the reference genome data.
   params {
     genomes {
       't7' {
         bowtie2_index     = "<path_to_refgenie_genomes>/alias/t7/bowtie2_index/default/t7"
         fasta             = "<path_to_refgenie_genomes>/alias/t7/fasta/default/t7.fa"
       }
     }
   }
   ```

4. Run your pipeline with the required genome. For example:

   :::bash
   nextflow run nf-core/<pipeline_name> --genome t7
   :::

### Handle asset name mismatches

When Refgenie asset names differ from nf-core expectations:

1. Create `alias_translations.yaml` in your Refgenie config directory.
2. Map Refgenie aliases to nf-core parameter names. For example:

   ```yaml title="alias_translations.yaml"
   ensembl_gtf: gtf
   star_index: star
   ```

3. Save the file.

:::note
Configuration automatically uses the translations.
:::
