---
title: Run your first pipeline
subtitle: Learn how to run the nf-core/demo pipeline
shortTitle: Run your first pipeline
weight: 3
---

With Nextflow and your software dependency manager installed, you are ready to run your first nf-core pipeline!
Running a test pipeline is ideal for verifying your environment is correctly configured, learning how nf-core pipelines work, and testing your installation.

## nf-core/demo

[`nf-core/demo`](https://nf-co.re/demo) is a lightweight demonstration pipeline designed to introduce running nf-core pipelines. It performs basic quality control and processing on genomic sequencing data, showcasing the structure and best practices followed by all nf-core pipelines.

The `nf-core/demo` pipeline runs three processes:

1. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
   - Quality control checks on raw sequencing reads
2. Adapter and quality trimming ([seqtk trim](https://github.com/lh3/seqtk))
   - Removes adapter sequences and low-quality bases
3. QC report aggregation ([MultiQC](https://seqera.io/multiqc/))
   - Combines QC metrics into a single interactive report

It is purposefully lightweight and should run within a few minutes on most systems.

## Quick start

The following sections cover common operations for getting started with `nf-core/demo`.

### Use test data

To run `nf-core/demo` with test data:

1. Run `nf-core/demo` with the built-in test profile and a profile to manage software dependencies:

   ```bash
   nextflow run nf-core/demo -profile test,<docker> --outdir results
   ```

   Replace `<docker>` with your preferred software dependency and compute environment manager.

   :::tip
   Common configure software dependencies and compute environments include:
   - `docker`: Uses Docker containers
   - `singularity`: Uses Singularity containers
   - `conda`: Uses Conda environments

   Combine profiles with commas: `-profile test,docker`.
   For a full list of supported software management tools, see [Software dependencies](./environment_setup/software-dependencies.md).
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

### Use your own data

To run `nf-core/demo` with your own sequencing data:

1. Create a comma-separated table called `samplesheet.csv` with the following columns:

   | Column    | Description                                                         |
   | --------- | ------------------------------------------------------------------- |
   | `sample`  | Unique sample identifier                                            |
   | `fastq_1` | Path to gzipped first-read FastQ file                               |
   | `fastq_2` | Path to gzipped second-read FastQ file (leave empty for single-end) |

   Example:

   ```csv
   sample,fastq_1,fastq_2
   SAMPLE1,</path/to/sample1_R1.fastq.gz>,</path/to/sample1_R2.fastq.gz>
   SAMPLE2,</path/to/sample2_R1.fastq.gz>,
   ```

1. Run `nf-core/demo` with your `samplesheet.csv` as input:

   ```bash
   nextflow run nf-core/demo \
     -profile <docker> \
     --input samplesheet.csv \
     --outdir results
   ```

   Replace `<docker>` with your preferred software dependency and compute environment manager.

### Use parameter files

To run `nf-core/demo` with a Nextflow parameters file:

1. Create `params.json` with the following contents:

   ```json
   {
     "input": "samplesheet.csv",
     "outdir": "results"
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
   They can also be used to publish alongside scientific publications for reproducibility purposes - but make sure to exclude hardcoded file paths specific to your system!
   :::

### Resume interrupted runs

To resume an interrupted run (for example due to an error, or manual cancellation), add the `-resume` flag:

```bash
nextflow run nf-core/demo \
  -profile <docker> \
  --input samplesheet.csv \
  --outdir results \
  -resume
```

Replace `<docker>` with your preferred software dependency and compute environment manager.

Only steps with changed inputs will re-execute, saving time and computational resources.

:::tip
Try pressing <kbd>Ctrl</kbd>+<kbd>c</kbd> during midway through pipeline run to simulate an interruption.
Then re-run the same command with `-resume` to see how it picks up from where it left off, due to the presence of 'cached' pipeline steps.
:::

### Version control

To specify a pipeline version, add the `-r` option with a version, revision, or commit ID:

```bash
nextflow run nf-core/demo \
  -r 1.0.2 \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

:::tip
Version information is automatically recorded in pipeline reports ine `pipeline_info/` results subdirectory.
:::

## Next steps

Now that you've successfully run your first nf-core pipeline:

- Browse the [nf-core pipeline catalog](https://nf-co.re/pipelines) to find workflows for your research area
- Learn to [adjust resource requirements and parameters](../running/configuration/overview.md) for your infrastructure
- Join the [nf-core Slack](https://nf-co.re/join/slack) community
