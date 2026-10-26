<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-example_logo_dark.png">
    <img alt="nf-core/example" src="docs/images/nf-core-example_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/example/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/example/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/example/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/example/actions/workflows/linting.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

## :warning: Deprecated parameters

The `--old_param` option was removed in version 2.0.0. See the [changelog](CHANGELOG.md#v200) for details.

## Introduction

**nf-core/example** is a bioinformatics pipeline that takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and alignment, and produces an extensive QC report.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/example_metro_map_dark.png">
  <img alt="Metro map" src="docs/images/example_metro_map.png" width="800">
</picture>

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment ([`STAR`](https://github.com/alexdobin/STAR) or [`HISAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml))
4. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A3_S3_L002_R1_001.fastq.gz,
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/example \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](docs/usage.md) and the [parameter documentation](https://nf-co.re/example/parameters).

| Parameter         | Description                              | Default  |
| ----------------- | ---------------------------------------- | -------- |
| `--input`         | Path to comma-separated samplesheet file | -        |
| `--outdir`        | The output directory                     | -        |
| `--genome`        | Name of iGenomes reference               | `GRCh38` |
| `--aligner`       | Alignment tool to use                    | `star`   |
| `--skip_trim`     | Skip the adapter trimming step           | `false`  |
| `--multiqc_title` | MultiQC report title                     | -        |

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/example/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](docs/output.md#multiqc).

> [!TIP]
> Images in the output docs such as ![MultiQC](docs/images/mqc_fastqc_counts.png) are rewritten to raw GitHub URLs.

## Credits

nf-core/example was originally written by the nf-core community.

Thanks to the following people for their extensive assistance in the development of this pipeline: see [CONTRIBUTORS](.github/CONTRIBUTING.md).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#example` channel](https://nfcore.slack.com/channels/example) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
