---
title: Reference Genomes
subtitle: How reference genomes are handled in nf-core
---

Many nf-core pipelines need a reference genome for alignment, annotation, or similar purposes.

## Illumina AWS iGenomes

:::warning
The transcriptome and GTF files in iGenomes are vastly out of date with respect to current annotations from Ensembl e.g. human iGenomes annotations are from Ensembl release 75, while the current Ensembl release is 108. Please consider downloading and using a more updated version of your reference genome as outlined in the next section.
:::

:::warning
The GRCh38 iGenomes assembly is from the NCBI and not Ensembl and as such there are some discrepancies in the way that the annotation is defined that may cause problems when running certain pipelines e.g. [nf-core/rnaseq#460](https://github.com/nf-core/rnaseq/issues/460). If you would like to use the latest soft-masked Ensembl assembly for GRCh38 instead please see the next section.
:::

To make the use of reference genomes easier, Illumina has developed a centralised resource called [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).
The most commonly used reference genome files are organised in a consistent structure for multiple genomes.

We have uploaded a copy of iGenomes onto AWS S3 and nf-core pipelines are configured to use this by default.
AWS iGenomes is hosted by Amazon as part of the [Registry of Open Data](https://registry.opendata.aws/aws-igenomes/) and are free to use. For more information about the AWS iGenomes, see [https://ewels.github.io/AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/).

All AWS iGenomes paths are specified in pipelines that support them in [`conf/igenomes.config`](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/igenomes.config#L14-L26). By default, the pipeline will automatically download the required reference files when you run the pipeline and supply an appropriate genome key (eg. `--genome GRCh37`). The pipeline will only download what it requires e.g. the nf-core/rnaseq pipeline will download the `star` indices specified in `conf/igenomes.config` but not the `bismark` index because that is something specific to the nf-core/methylseq pipeline. Genome asset related parameters required by nf-core pipelines are typically defined in the [main script](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/main.nf#L20-L28) for DSL2 pipelines. When using AWS iGenomes, for convenience, when a reference asset is available for direct download, these parameters are essentially auto-populated based on what is defined in `conf/igenomes.config` when you provide the `--genome` parameter. Downloading reference genome files takes time and bandwidth so, if possible, we recommend storing a local copy of your relevant iGenomes references which is outlined [here](https://ewels.github.io/AWS-iGenomes/).

To use a **local** version of iGenomes, the variable `params.igenomes_base` must be set to the path of the local iGenomes folder to reflect what is defined in [`conf/igenomes.config`](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/igenomes.config#L14-L26). Additional information on how to set the local `--igenomes_base` parameter can be found [here](/docs/usage/troubleshooting#using-a-local-version-of-igenomes).

To get the version of the annotation used by AWS iGenomes you can download the `README` file specified in the [`conf/igenomes.config`](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/igenomes.config#L22) via the AWS CLI e.g.

```bash
$ aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt .
download: s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt to ./README.txt

$ cat README.txt
The contents of the annotation directories were downloaded from Ensembl on: July 17, 2015.

Gene annotation files were downloaded from Ensembl release 75. SmallRNA annotation files were downloaded from miRBase release 21.
```

## Custom genomes

As mentioned in the section above, most of the required genome assets will be defined in the [main script](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/main.nf#L20-L28) for DSL2 nf-core pipelines. If you are unable to use the AWS iGenomes references, you can still supply reference genome parameters on the command line or via a `-params-file` in `yaml` or `json` format.

Using GRCh38 as an example, we can download the latest FASTA and GTF files using a simple bash script like below:

```bash
#!/bin/bash

VERSION=108
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz
```

Most genomics nf-core pipelines are able to start from just a FASTA and GTF file and create any downstream reference assets as part of the pipeline execution e.g. genome indices, intervals files etc. To avoid having to recreate these assets every time you run the pipeline you can use the `--save_reference` parameter that will save the indices, interval files etc in the results directory for you to move and store in a more central location for re-use with future pipeline runs. Using nf-core/rnaseq as an example [see docs](https://nf-co.re/rnaseq/output#reference-genome-files):

```bash
nextflow run \
    nf-core/rnaseq \
    --input samplesheet.csv \
    --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
    --gtf Homo_sapiens.GRCh38.108.gtf.gz \
    --save_reference \
    <OTHER_PARAMETERS>
```

Any downstream reference assets will be published in the results folder. For example, if you ran the nf-core/rnaseq pipeline in the step above with default options then the STAR index will be created and stored in the `<RESULTS_DIR>/genome/index/star` folder. Once you have moved the reference files to a central location so they are persistently available you can remove the `--save_reference` parameter and now explicitly override the relevant parameters via the CLI or a `-params-file` in `yaml` or `json` format. This will save having to re-create the genome indices and other assets over and over again which will be cost and time expensive.

```bash
nextflow run \
    nf-core/rnaseq \
    --input samplesheet.csv \
    --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
    --gtf Homo_sapiens.GRCh38.108.gtf.gz \
    --star_index /path/to/moved/star/directory/ \
    --gene_bed /path/to/moved/genes.bed \
    <OTHER_PARAMETERS>
```

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

:::note
You can also use [custom assets](http://refgenie.databio.org/en/latest/custom_assets/).
:::

3. Run your pipeline, specifying the required genome.

```bash
nextflow run nf-core/<PIPELINENAME> --genome t7  # ..rest of normal pipeline flags..
```

Please refer to [Refgenie documentation](http://refgenie.databio.org/en/latest/) for further information.

### How to handle Refgenie assets having different aliases than nf-core

A Refgenie server contains assets with established aliases, which can differ from the ones required by an nf-core pipeline.
For example, the asset for an ensemble index on the default Refgenie server is called [ensembl_gtf](http://refgenomes.databio.org/v3/assets/splash/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/ensembl_gtf?tag=default), while the same asset is called [gtf](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/igenomes.config#L33) in nf-core pipelines.

To address this, you can create a file with translations that will be used to generate the genomes configuration file with the appropriate names.

1. When you init Refgenie, you provide a path to a `genomes_config.yml` file with the argument `-c` or setting the environment variable `$REFGENIE`.
   In the same directory, create a file called `alias_translations.yaml`.

2. `alias_translations.yaml` must contain the equivalences of asset aliases in yaml format.
   Keys correspond to the name of refgenie server aliases while values correspond to the name of the respective nf-core pipeline aliases.
   For example:

   ```
   ensembl_gtf: gtf
   star_index: star
   ```

3. Pull your assets as usual. The asset aliases will be translated automatically in your `~/.nextflow/nf-core/refgenie_genomes.config` file.
