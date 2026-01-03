---
title: Testing
subtitle: Testing specifications for nf-core Nextflow DSL2 subworkflows
weight: 6
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Scope of testing

Tests for subworkflows SHOULD be designed to be executable within the nf-core/modules GitHub repository CI with example test data.

Tests for subworkflows MUST, at a minimum, run on the GitHub repository CI with a stub test that replicates the generation of (empty) output files.

Subworkflows tests do not necessarily need to execute 'standalone' (that is, run outside the nf-core/modules repository). For example, they do not need to be executable within a pipeline repository.

:::info{title="Rationale" collapse}
Some modules may require upstream modules or subworkflows to generate input files for the new module under construction if it is not possible or reasonable to upload those test data files to nf-core/test-datasets.

If the test was to work 'standalone,' the pipeline would need to include all these upstream modules/subworkflows just to execute the module test—even if those modules are not used within the pipeline itself. This would lead to a lot of file 'pollution' within the pipeline repository.

Subworkflows installed in the pipeline should already be tested to work correctly within the context of the pipeline with workflow- or pipeline-level tests. Thus, it is considered unnecessary to duplicate subworkflow tests again.
:::

:::note
CI tests for nf-core modules, subworkflows, or pipeline are **not** required to produce _meaningful_ output.

The main goal for nf-core CI tests are to ensure a given tool 'happily' executes without errors.

It is OK for a test to produce nonsense output, or find 'nothing', as long as the tool does not crash or produce an error.
:::

## All output channels must be tested

All output channels SHOULD be included in the nf-test snapshot file, or at a minimum, the files MUST be verified to exist.

## Tags

Tags MUST be specified for any dependent modules to ensure changes to upstream modules re-trigger tests for the subworkflow.

```groovy
tag "subworkflows"
tag "subworkflows_nfcore"
tag "<subworkflow_name>"
tag "<tool>" // Add each tool as a separate tag
tag "<tool>/<subtool>" // Add each subtool as a separate tag
```

## `assertAll()`

The `assertAll()` function MUST be used to specify an assertion, and at least one success assertion and versions MUST be included in the snapshot.

## Assert each type of input and output

A test and assertions SHOULD be written for each type of input and output.

Use [different assertion types](/docs/contributing/nf-test/assertions) if a straightforward `workflow.out` snapshot is not feasible.

:::tip
Always check the snapshot to ensure that all outputs are correct!
For exmaple, make sure there are no md5sums representing empty files.
:::

## Test names

Test names SHOULD describe the test dataset and configuration. For example:

```groovy
test("homo_sapiens - [fastq1, fastq2] - bam")
test("sarscov2 - [ cram, crai ] - fasta - fai")
test("Should search for zipped protein hits against a DIAMOND db and return a tab separated output file of hits")
```

## Input data

Input data SHOULD be referenced with the `modules_testdata_base_path` parameter:

```groovy
file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true)
```

:::info
CI tests for nf-core modules, subworkflows, or pipeline are **not** required to produce _meaningful_ output.

The main goal for nf-core CI tests are to ensure a given tool 'happily' executes without errors.

It is OK for a test to produce nonsense output, or find 'nothing', as long as the tool does not crash or produce an error.

Existing test-data from the modules branch of [nf-core/test-datasets](https://github.com/nf-core/test-datasets) SHOULD be reused as far as possible to reduce the size of our test dataset repository.

New test data SHOULD only be uploaded to nf-core/test-datasets if there is absolutely no other option within the existing test-data archive.
:::

## Configuration

Subworkflow nf-tests SHOULD use a single `nextflow.config` to supply `ext.args` to a subworkflow. Define them in the `when` block of a test under the `params` scope.

```groovy {4-7} title="main.nf.test"
config './nextflow.config'

when {
  params {
    moduleA_args = '--extra_opt1 --extra_opt2'
    moduleB_args = '--extra_optX'
  }
  process {
    """
    input[0] = [
      [ id:'test1', single_end:false ], // meta map
      file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)
    ]
    """
  }
}
```

```groovy {3,6} title="nextflow.config"
process {
  withName: 'MODULEA' {
    ext.args = params.moduleA_args
  }
  withName: 'MODULEB' {
    ext.args = params.moduleB_args
  }
}
```

Do not add other settings to this file.

:::tip
Supply the config only to the tests that use `params`, otherwise define `params` for every test including the stub test.
:::

## Skipping CI test profiles

If a subworkflow does not support a particular test profile, skip it by adding the path to the corresponding section in `.github/skip_nf_test.json`.

:::note
Keep the file sorted alphabetically.
:::
