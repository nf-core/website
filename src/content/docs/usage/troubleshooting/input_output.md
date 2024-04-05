---
title: Input and output errors
subtitle: How to troubleshoot common mistakes and issues
---

## Input and Output errors

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

### Direct input

:::warning
This section mostly refers to DSL1 pipelines! Most DSL2 pipelines now use samplesheet inputs rather than direct read inputs.
:::

Or when you're using a input method like `--input '/<path>/<to>/*_fq.gz'`, but only pick up one file, or only one file per pair being processed during the run, please note the following:

1. [The path must be enclosed in quotes (`'` or `"`)](#output-for-only-a-single-sample-although-i-specified-multiple-with-wildcards)
2. The path must have at least one `*` wildcard character i.e. following a ['glob' pattern](https://en.wikipedia.org/wiki/Glob_%28programming%29). This is even if you are only running one paired end sample.
   - A description of valid pattern matching can be seen [here](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) for java and [here](https://www.nextflow.io/docs/latest/channel.html?highlight=glob#frompath) for Nextflow
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
   - This notation is interpreted by Nextflow to mean anything with the same string other than R1 and R2 in the file name, will be be assumed to be a pair of files.
4. If you are running single-end data make sure to specify `--single_end`
5. [Your data should be organised in a 'tidy' manner](#data-organization)

A few examples are as follows:

- Running with a single, single-end FASTQ file as input (this will produce output files for this sample only)

  ```bash
  nextflow run nf-core/<pipeline> --input 'my_data.fastq.gz` --single_end
  ```

- Running multiple single-end FASTQ files as input using a wildcard glob pattern. This will find all files in the directory beginning with `my_`, and ending in `.fastq.gz`, with each file with any other characters between those two strings being considered distinct samples (and will produce output files for each of the multiple input files).

  ```bash
  nextflow run nf-core/<pipeline> --input 'my_*.fastq.gz` --single_end
  ```

- Running multiple paired-end FASTQ files as input using wildcard and grouping glob patterns. This will find all files in the directory beginning with `my_`, and ending in `.fastq.gz`, with each file with any other characters between those two strings being considered distinct samples. However, any pair of file names that are exactly the same other than `R1` and `R2` will be grouped together, and processed as related files. You will in most cases get output files for each distinct file, but with the `*{R1,R2}` syntax, R1 and R2 pairs are collapsed into one.

  ```bash
  nextflow run nf-core/<pipeline> --input 'my_*{R1,R2}.fastq.gz`
  ```

Note that if your sample name is "messy" then you have to be very particular with your glob specification (see point 2 above). A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` will not give you what you want, whilst `*{R1,R2}*.gz` will.

Please also note that genomics pipelines with 'direct input' can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. Furthermore, you have paired-end sequencing data mode please make sure that your files are sensibly named so that they can be properly paired e.g. (`--input '*_{R1,R2}.fastq.gz'`).

### Output for only a single sample although I specified multiple with wildcards

You must specify paths to files in quotes, otherwise your _shell_ (e.g. bash) will evaluate any wildcards (\*) rather than Nextflow.

For example:

```bash
nextflow run nf-core/<pipeline> --input /path/to/sample_*/*.fq.gz
```

Maybe evaluated by your shell as:

```bash
nextflow run nf-core/<pipeline> --input /path/to/sample_1/sample_1.fq.gz /path/to/sample_1/sample_1.fq.gz /path/to/sample_1/sample_1.fq.gz
```

And Nextflow will only take the first path after `--input`, ignoring the others.

On the other hand, encapsulating the path in quotes will allow _Nextflow_ to evaluate the paths.

```bash
nextflow run nf-core/<pipeline> --input "/path/to/sample_*/*.fq.gz"
```

### Sample sheet input

If you are using a sample sheet or TSV input method, check there is not a mistake or typo in the path in a given column. Common mistakes are a trailing space at the end of the path, which can cause problems.

### Using a local version of iGenomes

The iGenomes config file uses `params.igenomes_base` to make it possible to use a local copy of iGenomes. However, as custom config files are loaded after `nextflow.config` and the `igenomes.config` has already been imported and parsed, setting `params.igenomes_base` in a custom config file has no effect and the pipeline will use the s3 locations by default. To overcome this you can specify a local iGenomes path by either:

- Specifying an `--igenomes_base` path in your execution command.

```bash
nextflow run nf-core/<pipeline> --input <input> -c <config> -profile <profile> --igenomes_base <path>/<to>/<data>/igenomes
```

- Specifying the `igenomes_base` parameter in a parameters file provided with `-params-file` in `yaml` or `json` format.

```bash
nextflow run nf-core/<pipeline> -profile <profile> -params-file params.yml
```

Where the `params.yml` file contains the pipeline params:

```yaml title="params.yml"
input: '/<path>/<to>/<data>/input'
igenomes_base: '/<path>/<to>/<data>/igenomes'
```
