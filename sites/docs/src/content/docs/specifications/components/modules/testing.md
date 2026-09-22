---
title: Testing
subtitle: Test modules
markdownPlugin: addNumbersToHeadings
shortTitle: Testing
weight: 8
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

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

## GPU tests

Modules that support both CPU and GPU modes SHOULD include a separate GPU test file (`main.gpu.nf.test`).
GPU-only modules MAY use a single test file.

GPU tests MUST be tagged with `"gpu"` or `"gpu_highmem"` so the GPU CI workflow discovers and runs them on GPU-enabled runners.
GPU tests SHOULD include a `nextflow.gpu.config` that sets `accelerator = 1` on the process.
GPU tests SHOULD use the same assertions as the CPU tests to verify that GPU and CPU modes produce equivalent results, and SHOULD include both a real test and a stub test.

For runner instance types, the GPU concurrency caveat under Singularity, and a worked example, see the [GPU-capable modules](/docs/developing/components/gpu-modules) guide.

## Stub tests

A stub test MUST be included for the module.

## Tags

Tags MUST be specified for any dependent modules to ensure changes to upstream modules will re-trigger tests for the current module.

```groovy
tag "modules"
tag "modules_nfcore"
tag "<tool>"
tag "<tool>/<subtool>" // Only if there is a subtool
tag "<dependent-tool>/<dependent-subtool>" // Only if there is a tool this module depends on
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
Give `module_args` a default in the config's `params` scope, and apply the config to each test individually rather than once at the top of the file, so a test can override the default in its own `params` block when it needs to:

```groovy {2-4} title="nextflow.config"
params {
  module_args = ''
}

process {
  withName: 'MODULE' {
    ext.args = params.module_args
  }
}
```

No other settings should go into this file.

```groovy {2,6-8} title="main.nf.test"
test("my_tool - custom args") {
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
}

test("my_tool - stub") {
  config './nextflow.config'
  options '-stub'
  when {
    process {
      """
      input[0] = [
        [ id:'test1', single_end:false ], // meta map
        file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)
      ]
      """
    }
  }
}
```

Because `module_args` already has a default in the config, tests that don't need custom args (like the stub test above) don't need to repeat `module_args = ''` — they inherit the default as soon as they apply `config './nextflow.config'`.

:::info
Modules in pipelines are frequently configured with dynamic inputs. Test parameters do not support this. For example,

```groovy {3-4} title="nextflow.config"
process {
  withName: 'MODULE' {
    ext.args = { "--sample ${meta.id}" }
    ext.prefix = { "${meta.id}_prefix" }
  }
}
```

would be implemented as follows:

```groovy {2,6-8} title="main.nf.test"
test("my_tool - dynamic args") {
  config './nextflow.config'
  when {
    params {
      module_args = '--sample test1' // `meta.id` is replaced with the value it would take once dynamically resolved
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
}
```

```groovy {2-6} title="nextflow.config"
params {
  module_args = ''
}
process {
  withName: 'MODULE' {
    ext.args = params.module_args
    ext.prefix = { "${meta.id}_prefix" } // Dynamic prefix configuration remains in the config
  }
}
```

:::

## Skipping CI test profiles

If a module does not support a particular test profile, you can skip it by adding the path to corresponding section in `.github/skip_nf_test.json`.

:::note
Sort the file alphabetically.
:::
