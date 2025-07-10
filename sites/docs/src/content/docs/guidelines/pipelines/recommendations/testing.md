---
title: Testing
subtitle: Use nf-test to validate pipeline
weight: 210
---

Pipelines should use nf-test for valid pipeline testing.

All pipelines should support nf-test and use it for pipeline level tests.
Additional tests for local subworkflows and local modules are recommended.
Modules and subworkflows from nf-core/modules should already include nf-tests which can also be used within the pipeline.

Pipeline level tests can facilitate more reliable and reproducible pipelines by ensuring the pipeline produces identical results with every run.
You must add pipeline tests that work with `-profile test` and you should reuse this profile within one nf-test.

### Pipeline nf-test overview and structure

Within the base directory of the repository, there is a configuration file for nf-test, named `nf-test.config`.
This will set the following options:

- Set the `testsDir` to the base of the repository so it includes all files
- Set the `workDir` to `NFT_WORKDIR` if it is set, otherwise it will default to `.nf-test`
- Add an additional configuration file specific for nf-test located in `tests/nextflow.config`
- Set the default profile(s) for nf-test to include `test` (this can be overridden on the command line or in the test files `*.nf.test`. See nf-test docs about [managing profiles](https://www.nf-test.com/docs/configuration/#managing-profiles).)
- Set filenames that should trigger a full test run when modified (`nextflow.config`, `nf-test.config`, `conf/test.config`, `tests/nextflow.config`)
- Load the necessary plugins (`nft-utils` is used for parsing the output of the pipeline, and other plugins can be added as needed ie `nft-bam`, `nft-vcf`...)

The `tests/` folder contains all nf-test related files.

The `.nftignore` file that lists files that should be ignored by nf-test using the `nft-utils` plugin.
Within the nf-test specific configuration file (`tests/nextflow.config`), you can add specific requirements for running nf-tests but this should not include parameters or options as these should be available in all contexts.
The default pipeline level test (`main.nf.test`) is used to test the pipeline with the profile test, and should be run once to capture the output of the pipeline.
That will generate snapshots (`main.nf.test.snap`) that will be used to compare the output of the pipeline with the expected output.
Additional nf-test file must contain related pipeline tests which tests for a single set up.
It makes it easier to launch related pipeline tests at once.
Each nf-test file should be named after the set-up it tests in the following format:

```tree
tests
├─ alignment-pair-end.nf.test
├─ alignment-single-end.nf.test
└─ default.nf.test
```

### Pipeline nf-tests additional guidance

The same guidelines for test profiles, test data and nf-test also apply to pipeline tests.
In addition, the following guidelines apply:

- To ensure all output files are caught, the `params.outdir` should be set to the nf-test variable `outputDir`
- To ensure all output files are caught, it is highly recommended to use the `nft-utils` plugin to parse the output of the pipeline.
- To ensure a proper execution of the pipeline is captured, it is recommended to capture the number of successful tasks.
- To ensure that all output files are captured, it is recommended to run the tests once, and rerun to verify they're stable.
  - If not feasible, such files can be added to the `.nftignore` file.

```groovy title="main.nf.test"
nextflow_pipeline {

    name "Test pipeline"
    script "../main.nf"
    tag "pipeline"

    test("-profile test") {

        when {
            params {
                outdir = "$outputDir"
            }
        }

        then {
            // stable_name: All files + folders in ${params.outdir}/ with a stable name
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
            // stable_path: All files in ${params.outdir}/ with stable content
            def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
            assertAll(
                { assert workflow.success},
                { assert snapshot(
                    // Number of successful tasks
                    workflow.trace.succeeded().size(),
                    // pipeline versions.yml file for multiqc from which Nextflow version is removed because we tests pipelines on multiple Nextflow versions
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_{{ short_name }}_software_mqc_versions.yml"),
                    // All stable path name, with a relative path
                    stable_name,
                    // All files with stable contents
                    stable_path
                ).match() }
            )
        }
    }
}
```

`stable_name` capture all of the stable files and folders in the output directory.
`stable_path` capture all of the stable files in the output directory, ignoring the files listed in the `.nftignore` file.

Additional plugins can be used if files with unstable content can be read using a nf-test plugin.
For example for BAM and CRAM files the `nft-bam` plugin is used, while the `nft-vcf` plugin is also used, both in conjonction with `nft-utils`:

```groovy
nextflow_pipeline {

    name "Test pipeline"
    script "../main.nf"
    tag "pipeline"

    test("-profile test") {

        when {
            params {
                outdir = "$outputDir"
            }
        }

        then {
            // stable_name: All files + folders in ${params.outdir}/ with a stable name
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
            // stable_path: All files in ${params.outdir}/ with stable content
            def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
            // bam_files: All bam files
            def bam_files  = getAllFilesFromDir(params.outdir, include: ['**/*.bam'])
            // cram_files: All cram files
            def cram_files  = getAllFilesFromDir(params.outdir, include: ['**/*.cram'])
            // Fasta file for cram verification with nft-bam
            def fasta_base  = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
            def fasta       = fasta_base + 'genomics/homo_sapiens/genome/genome.fasta'
            // vcf_files: All vcf files
            def vcf_files  = getAllFilesFromDir(params.outdir, include: ['**/*.vcf.gz'])
            assertAll(
                { assert workflow.success},
                { assert snapshot(
                    // Number of successful tasks
                    workflow.trace.succeeded().size(),
                    // pipeline versions.yml file for multiqc from which Nextflow version is removed because we tests pipelines on multiple Nextflow versions
                    removeNextflowVersion("$outputDir/pipeline_info/nf_core_{{ short_name }}_software_mqc_versions.yml"),
                    // All stable path name, with a relative path
                    stable_name,
                    // All files with stable contents
                    stable_path,
                    // All bam files
                    bam_files.collect { file -> [file.getName(), bam(file.toString()).readsMD5] },
                    // All cram files
                    cram_files.collect { file -> [file.getName(), cram(file.toString(), fasta).readsMD5] },
                    // All vcf files
                    vcf_files.collect { file -> [file.getName(), path(file.toString()).vcf.variantsMD5] }
                ).match() }
            )
        }
    }
}
```

The lists generated by nft-utils can be empty, and are therefore including an empty line with trailing whitespace in the snapshot.
It should be removed from the test, but a ternary operator can be used to avoid this and display a string instead of an empty list:

```groovy
stable_path.isEmpty() ? 'No files with stable names' : stable_path
```

This could be useful when the same structure for the test is copied over from one test to another.
