---
title: Run your first pipeline
subtitle: Run the nf-core/demo pipeline
shortTitle: Your first pipeline
weight: 3
---

With Nextflow and your software dependency manager installed, you are ready to run your first nf-core pipeline.

## nf-core/demo

[`nf-core/demo`](https://nf-co.re/demo) is a lightweight demonstration pipeline designed to introduce running nf-core pipelines. It performs basic quality control and processing on sequencing data, showcasing the structure and best practices followed by all nf-core pipelines.

It is ideal for:
- Verifying your environment is correctly configured
- Learning how nf-core pipelines work
- Testing your installation

The pipeline runs three processes:

1. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    - Quality control checks on raw sequencing reads
2. Adapter and quality trimming ([seqtk trim](https://github.com/lh3/seqtk))
    - Removes adapter sequences and low-quality bases
3. QC report aggregation ([MultiQC](https://seqera.io/multiqc/))
    - Combines QC metrics into a single interactive report

## Quick start with test data

To get started quickly:

1. Run `nf-core/demo` with the built-in test profile and a profile to manage software dependencies:

    ```bash
    nextflow run nf-core/demo -profile test,docker --outdir results
    ```

    :::tip
    Common configure software dependencies and compute environments include:

    - `docker`: Uses Docker containers
    - `singularity`: Uses Singularity containers
    - `conda`: Uses Conda environments

    Combine profiles with commas: `-profile test,docker`
    :::

1. Once complete, view your `results/` directory:

    ```console
    results/
    ├── fastqc/          # Raw FastQC reports for each sample
    ├── fq/              # Trimmed FastQ files
    ├── multiqc/         # Aggregated QC report
    │   └── multiqc_report.html
    └── pipeline_info/   # Execution logs and metadata
    ```

## Running with your own data

To run `nf-core/demo` with your own sequencing data:

1. Create `samplesheet.csv` with the following contents:

    | Column | Description |
    |--------|-------------|
    | `sample` | Unique sample identifier |
    | `fastq_1` | Path to gzipped first-read FastQ file |
    | `fastq_2` | Path to gzipped second-read FastQ file (leave empty for single-end) |

    Example:

    ```csv
    sample,fastq_1,fastq_2
    SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
    SAMPLE2,/path/to/sample2_R1.fastq.gz,
    ```

1. Run `nf-core/demo` with your `samplesheet.csv` as input:

    ```bash
    nextflow run nf-core/demo \
      -profile <docker> \
      --input samplesheet.csv \
      --outdir results
    ```

    Replace `<docker>` with your preferred software dependency and compute environment manager.


## Using parameter files

To run `nf-core/demo` with a parameters file:

1. Create `params.json` with the following contents:

    ```json
    params {
      input = 'samplesheet.csv' // This requires the samplesheet.csv created above
      outdir = 'results'
    }
    ```

1. Run `nf-core/demo` with you `params.json` file:

    ```bash
    nextflow run nf-core/demo \
      -profile <docker> \
      -params-file params.json
    ```

    Replace `<docker>` with your preferred software dependency and compute environment manager.

    :::tip
    Parameter files make it easier to rerun pipelines with the same settings and keep track of your analysis configurations.
    :::

## Resuming interrupted runs

Nextflow caches completed steps. If a run is interrupted, use `-resume` to continue:

```bash
nextflow run nf-core/demo \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

Only steps with changed inputs will re-execute, saving time and computational resources.

## Version control for reproducibility

For reproducible analyses, specify a pipeline version with `-r`:

```bash
nextflow run nf-core/demo \
  -r 1.0.2 \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

:::tip
Version information is automatically recorded in pipeline reports. Check [available releases](https://github.com/nf-core/demo/releases) on GitHub.
:::

## Next steps

Now that you've successfully run your first nf-core pipeline:

- **Explore other pipelines**: Browse the [nf-core pipeline catalog](https://nf-co.re/pipelines) to find workflows for your research area
- **Customize configurations**: Learn to adjust resource requirements and parameters for your infrastructure
- **Join the community**: Get help and share experiences on [nf-core Slack](https://nf-co.re/join/slack)
