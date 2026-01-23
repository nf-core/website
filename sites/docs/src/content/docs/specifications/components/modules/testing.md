---
title: Testing
subtitle: Testing requirements and guidelines for nf-core modules
markdownPlugin: addNumbersToHeadings
shortTitle: Testing
weight: 8
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Scope of testing

Tests for modules SHOULD be executable within the nf-core/modules GitHub repository CI with example test data.

Tests for modules MUST, at a minimum, run on the GitHub repository CI with a stub test that replicates the generation of (empty) output files and a `versions` file.

Module tests do not necessarily need to be able to execute 'standalone', i.e., run outside the nf-core/modules repository.
For example, they don't need to be executable within a pipeline repository.

:::info{title="Rationale" collapse}
Some modules may require upstream modules to generate input files for the new module under construction if it is not possible or reasonable to upload those test data files to nf-core/test-datasets.

If the test were to work 'standalone,' the pipeline would need to include all these upstream modules to execute the module test—even if those modules are not used within the pipeline itself.
This would lead to a lot of file 'pollution' within the pipeline repository.

Modules installed in the pipeline should already be tested to work correctly within the context of the pipeline with workflow- or pipeline-level tests.
It is considered unnecessary to duplicate module tests again.
:::

:::note
CI tests for nf-core modules, subworkflows, or pipeline are **not** required to produce _meaningful_ output.

The main goal for nf-core CI tests are to ensure a given tool 'happily' executes without errors.

It is OK for a test to produce nonsense output, or find 'nothing', as long as the tool does not crash or produce an error.
:::

## Snapshots

Use only one snapshot per module test, which SHOULD contain all assertions present in this test.
Having multiple snapshots per test will make the snapshot file less readable.

The snapshot SHOULD include all output channels for each test, or at a minimum, it MUST contain some verification that the file exists.

By default, the `then` block of a test should contain this:

```groovy
assert snapshot(process.out).match()
```

When the snapshot is unstable, use another way to test the output files.
See [nf-test assertions](/docs/contributing/nf-test/assertions) for examples on how to do this.

## Stub tests

A stub test MUST be included for the module.

## Tags

Tags MUST be specified for any dependent modules to ensure changes to upstream modules will re-trigger tests for the current module.

```groovy
tag "modules"
tag "modules_nfcore"
tag "<tool>"
tag "<tool>/<subtool>" // Only if there is a subtool
tag "<dependent_tool>/<dependent_subtool>" // Only if there is a tool this module depends on
```

## `assertAll()`

Use the `assertAll()` function to specify an assertion, and there MUST be a minimum of one success assertion and versions in the snapshot.

## Assert each type of input and output

Include a test and assertions for each type of input and output.

Use [different assertion types](/docs/contributing/nf-test/assertions) if a straightforward `process.out` snapshot is not feasible.

:::tip
Always check the snapshot to ensure that all outputs are correct!
For example, make sure there are no md5sums representing empty files (with the exception of stub tests!).
:::

## Test names

Test names SHOULD describe the test dataset and configuration used. some examples below:

```groovy
test("homo_sapiens - [fastq1, fastq2] - bam")
test("sarscov2 - [ cram, crai ] - fasta - fai")
test("Should search for zipped protein hits against a DIAMOND db and return a tab separated output file of hits")
```

## Input data

Reference input data with the `modules_testdata_base_path` parameter:

```groovy
file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true)
```

:::info
CI tests for nf-core modules, subworkflows, or pipeline are **not** required to produce _meaningful_ output.

The main goal for nf-core CI tests are to ensure a given tool 'happily' executes without errors.

It is OK for a test to produce nonsense output, or find 'nothing', as long as the tool does not crash or produce an error.

Therefore, reuse existing test-data from the modules branch of [nf-core/test-datasets](https://github.com/nf-core/test-datasets) as far as possible to reduce the size of our test dataset repository.

Upload new test data to nf-core/test-datasets only if there is absolutely no other option within the existing test-data archive.
:::

## Configuration of ext.args in tests

Module nf-tests SHOULD use a single `nextflow.config` to supply `ext.args` to a module.
Define them in the `when` block of a test under the `params` scope.

```groovy {4-6} title="main.nf.test"
config './nextflow.config'

when {
  params {
    module_args = '--extra_opt1 --extra_opt2'
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

```groovy {3} title="nextflow.config"
process {
  withName: 'MODULE' {
    ext.args = params.module_args
  }
}
```

No other settings should go into this file.

:::tip
Supply the config only to the tests that use `params`, otherwise define `params` for every test including the stub test.
:::

## Skipping CI test profiles

If a module does not support a particular test profile, you can skip it by adding the path to corresponding section in `.github/skip_nf_test.json`.

:::note
Sort the file alphabetically.
:::
