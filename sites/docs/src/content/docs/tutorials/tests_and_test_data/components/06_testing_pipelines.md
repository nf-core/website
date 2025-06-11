---
title: "6. Testing Pipelines"
subtitle: End-to-end pipeline testing
weight: 60
---

# Pipeline Testing

Pipeline-level testing ensures that your entire nf-core pipeline works correctly from start to finish. 
As of nf-core/tools 3.3, pipeline-level nf-test tests have been added to the pipeline template to improve robustness and help developers catch issues early in the development process.

## Template Files

When you create a new nf-core pipeline or update an existing one, you'll find these new template files for pipeline testing:

```tree
nf-test.config       # The pipeline-level nf-test configuration file
tests/
├── .nftignore       # ignored files for nf-test
├── default.nf.test  # The default test for the pipeline, loading the parameter configuration from config/test.config
└── nextflow.config  # The nextflow configuration for the pipeline tests
```

## Pipeline Testing

nf-test allows you to test the complete pipeline end-to-end. The basic syntax for a pipeline test follows this structure:

```groovy
nextflow_pipeline {

    name "<NAME>"
    script "<PATH/TO/NEXTFLOW_SCRIPT.nf>"

    test("<TEST_NAME>") {

    }
}
```

The workflow object can be used in asserts to check its status, error messages or traces.

```groovy
// workflow status
assert workflow.success
assert workflow.failed
assert workflow.exitStatus == 0

// workflow error message
assert workflow.errorReport.contains("....")

// trace
//returns a list containing succeeded tasks
assert workflow.trace.succeeded().size() == 3

//returns a list containing failed tasks
assert workflow.trace.failed().size() == 0

//returns a list containing all tasks
assert workflow.trace.tasks().size() == 3
```

---

## Using nft-utils Plugin

The nft-utils plugin provides additional utilities for pipeline testing.

[nft-utils documentation](https://github.com/nf-core/nft-utils/blob/main/docs/usage.md)

### Key Concepts

#### Stable Files vs Stable Names

When testing pipelines, we need to distinguish between two types of outputs:

- **Stable Names**: Files and directories that have consistent names across different test runs, but whose content may vary (e.g., due to timestamps, random seeds, or system-specific information). These are tracked by their relative paths within the output directory.

- **Stable Files**: Files that have both consistent names AND consistent content across test runs. These files can be safely included in snapshots for content comparison.

#### .nftignore File

The `.nftignore` file works similarly to `.gitignore` but for nf-test snapshots. It allows you to exclude files with unstable content from being included in snapshot comparisons while still tracking their presence by name. This is essential for:

- Log files with timestamps
- HTML reports with dynamic content
- Binary files that may have slight variations
- Files containing system-specific paths or information

Files listed in `.nftignore` will be excluded from `stable_path` snapshots but can still be tracked in `stable_name` snapshots for structural verification.

### Installing and Configuring nft-utils

Add to your `nf-test.config`:

```groovy
config {
    plugins {
        load "nft-utils@0.0.4"
    }
}
```

### Essential Pipeline-Level Snapshot Pattern

This is the **recommended pattern** for comprehensive pipeline testing:

```groovy
nextflow_pipeline {
    name "Test Pipeline Default"
    script "main.nf"
    tag "pipeline"
    tag "default"

    test("Pipeline with stable snapshots") {
        when {
            params {
                input = 'tests/samplesheet.csv'
                outdir = "$outputDir"
                genome = 'GRCh37'
            }
        }
        then {
            // stable_name: All files + folders in ${params.outdir}/ with stable names
            def stable_name = getAllFilesFromDir(
                params.outdir,
                relative: true,
                includeDir: true,
                ignore: ['pipeline_info/*.{html,json,txt}', 'pipeline_info/execution_*.{html,txt}']
            )

            // stable_path: All files in ${params.outdir}/ with stable content
            def stable_path = getAllFilesFromDir(
                params.outdir,
                ignoreFile: 'tests/.nftignore'
            )

            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    // Number of successful tasks
                    workflow.trace.succeeded().size(),
                    // Pipeline versions.yml file with Nextflow version removed
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_*_software_mqc_versions.yml"),
                    // All stable file/folder names (relative paths)
                    stable_name,
                    // All files with stable contents
                    stable_path
                ).match() }
            )
        }
    }
}
```

### Setting up .nftignore for Unstable Files

Create `tests/.nftignore` to exclude files with unstable content but stable names. Here's an example:

```nftignore
.DS_Store
*/alignments/logs/*.txt
*/{alignments,deduplicated}/*.{bam,bam.bai}
*/deduplicated/logs/*.txt
*/{reports,summary}/*.{html,txt}
*/deduplicated/picard_metrics/*.txt
fastqc/*.html
fastqc/zips/*.zip
multiqc/*/multiqc_data/*.{log,json}
multiqc/*/multiqc_data/multiqc_fastqc.txt
multiqc/*/multiqc_data/multiqc_general_stats.txt
multiqc/*/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt
multiqc/*/multiqc_data/multiqc_software_versions.txt
multiqc/*/multiqc_data/multiqc_sources.txt
multiqc/*/multiqc_report.html
multiqc/*/multiqc_plots/{pdf,png,svg}/*.{pdf,png,svg}
pipeline_info/*.{html,json,txt,yml}
*/qualimap/bamqc/*/qualimapReport.html
*/qualimap/bamqc/**/*.{pdf,png,svg}
*/qualimap/bamqc/*/css/*
qualimap/**/images_qualimapReport/*
qualimap/**/raw_data_qualimapReport/*
trimgalore/fastqc/*.html
trimgalore/fastqc/zips/*.zip
trimgalore/logs/*.txt
enrichment_metrics/*
```

### Advanced nft-utils Usage

[Sarek is doing it](https://github.com/nf-core/sarek/blob/dev/tests%2Fdefault.nf.test#L38-L47)

```groovy
                { assert snapshot(
                    // Number of successful tasks
                    workflow.trace.succeeded().size(),
                    // pipeline versions.yml file for multiqc from which Nextflow version is removed because we tests pipelines on multiple Nextflow versions
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_sarek_software_mqc_versions.yml"),
                    // All stable path name, with a relative path
                    stable_name,
                    // All files with stable contents
                    stable_path.isEmpty() ? 'No stable content' : stable_path,
                    // All cram files
                    bam_files.isEmpty() ? 'No BAM files' : bam_files.collect { file -> file.getName() + ":md5," +  bam(file.toString()).readsMD5 },
                    // All cram files
                    cram_files.isEmpty() ? 'No CRAM files' : cram_files.collect { file -> file.getName() + ":md5," +  cram(file.toString(), fasta).readsMD5 },
                    // All vcf files
                    vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.getName() + ":md5," + path(file.toString()).vcf.variantsMD5 }
                ).match() }
```

### Best Practices for nft-utils

1. **Always set outdir to `$outputDir`** - ensures consistent test directory structure
2. **Use meaningful snapshot tags** - helps identify specific test scenarios
3. **Separate stable names from stable content** - provides better debugging when tests fail
4. **Include task count verification** - `workflow.trace.succeeded().size()` ensures execution completeness
5. **Remove version dependencies** - use `removeNextflowVersion()` for reproducible tests
6. **Maintain .nftignore files** - keep them updated as pipeline outputs evolve

### Troubleshooting nft-utils

#### When snapshots fail to match:

1. **Check .nftignore patterns** - ensure unstable files are properly excluded
2. **Verify wildcard patterns** - make sure they match your actual file structure
3. **Review relative vs absolute paths** - use `relative: true` for portable tests
4. **Update ignore patterns** - add new unstable files discovered during testing

## Next Steps

Continue to [nf-test Assertions](./07_assertions.md) to learn about comprehensive assertion patterns and verification techniques.
