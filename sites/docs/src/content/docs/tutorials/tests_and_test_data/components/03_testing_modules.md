---
title: "Testing nf-core/modules"
subtitle: Writing comprehensive tests for nf-core modules with nf-test
weight: 30
---

## Overview

The nf-test framework in the nf-core ecosystem enables comprehensive testing for processes, workflows, and pipelines.

This chapter covers the fundamentals of nf-core module testing, from basic syntax to advanced scenarios involving chained modules.

## Understanding snapshots

**What are snapshots?**

Snapshots are nf-test's way of capturing and validating the expected outputs of your tests. When you run a test for the first time, nf-test creates a `.snap` file containing checksums, file names, and other metadata that describes how a module's output should look like.

**How snapshot matching works:**

1. **First run**: nf-test generates a snapshot file (e.g., `main.nf.test.snap`) containing the expected outputs
2. **Subsequent runs**: nf-test compares current outputs against the stored snapshot
3. **Pass/Fail**: Test passes if outputs match the snapshot, fails if they differ

**Why snapshots are useful:**

- **Automated validation**: No need to manually check every output file
- **Regression detection**: Automatically catch when module behavior changes
- **Comprehensive checking**: Validates file content, structure, and metadata

## Your first simple test walkthrough with an example

Let's walk through creating a basic test step-by-step using the `cat` module as an example:

### Step 1: Examine the module

First, let's look at a simple example module (`modules/cat/main.nf`):

```groovy
process CAT {
    input:
    tuple val(meta), path(file_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions

    script:
    prefix = "${meta.id}.txt"
    """
    cat ${file_in} > ${prefix}
    """
}
```

This process simply concatenates input files and outputs a single file.

### Step 2: Write the test

Create `modules/cat/tests/main.nf.test`:

```groovy
nextflow_process {
    name "Test Process CAT"
    script "../main.nf"
    process "CAT"

    test("simple file concatenation") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test' ], // meta map with sample identifier
                    file('test_file.txt')  // input file
                ]
                """
            }
        }

        then {
            assert process.success               // Check process completed successfully
            assert snapshot(sanitizeOutput(process.out)).match() // Validate outputs match expected snapshot
        }
    }
}
```

**Understanding the test structure:**

The nf-test file is organized into several key scopes that define what and how to test:

- `nextflow_process { }` - **Test wrapper**: Declares this is a process-level test
- `name` - **Test description**: Human-readable name for the test suite
- `script` - **Module location**: Absolute/Relative path to the `main.nf` file containing the process. Relative paths should start from the location of the current test file.
- `process` - **Process name**: The specific process to test (must match the process name in `main.nf`)
- `test("...")` - **Individual test case**: Each test block defines one scenario to test
- `when { }` - **Test setup**: Where you define the input data for the process
- `then { }` - **Assertions**: Where you specify what and how to check that the process behaved as expected

**Inside the test components:**

Now that we understand the overall structure, let's look at the specific elements within the test:

**When section (input setup):**

- `[ id:'test' ]` - **meta map**: Contains sample metadata (ID, conditions, etc.)
- `file('test_file.txt')` - **input file**: The file to process

**Then section (assertions):**

- `process.success` - **success check**: Ensures the process completed without errors
- `snapshot(sanitizeOutput(process.out)).match()` - **output validation**: Compares all outputs to stored snapshot
- `sanitizeOutput()` - **snapshot sanitization**: Cleans process and workflow outputs by removing numbered keys, making snapshots more human-readable

### Step 3: Run the test

Once you've written the test, you can generate the initial snapshot by running the test using nf-core tools.

```bash
cd path/to/modules
nf-core modules test cat --profile docker
```

**What happens:**

1. nf-test creates a temporary work directory `.nf-test/`
2. Runs the CAT process with your test inputs
3. **First run**: Generates `tests/main.nf.test.snap` with output checksums
4. **Second run**: Validates outputs match the generated snapshot
5. Reports PASS/FAIL

:::info
Note that nf-core tools wraps around `nf-test` commands itself.
nf-core tools commands are not 1:1 equivalent with 'raw' nf-test commands.
:::

### Step 4: Examine the generated snapshot

After running, check `tests/main.nf.test.snap`:

```json
{
  "simple file concatenation": {
    "content": [
      {
        "0": [
          // First output channel (file_out)
          [{ "id": "test" }, "test.txt:md5,d41d8cd98f00b204e9800998ecf8427e"]
        ],
        "1": [
          // Second output channel (versions)
          "versions.yml:md5,c50aa59475ab20752531545a585f8f2d"
        ]
      }
    ]
  }
}
```

You can see the outputs of the module, with an associated `md5sum` which concisely describes the exact contents of the two output files: the concatenated file, and the nf-core `versions.yml` file.

Now you understand the full test cycle! The snapshot ensures your module produces consistent, expected outputs every time.

Any time you change the module - you can run the test again to check that you did not unexpectedly change an output that should not have been changed.

## Essential assertions

Tests use assertions to verify the expected output of the process specified in the `then` block.

You can specify multiple assertions to be evaluated together in a single test by specifying them within an `assertAll` block.

Nextflow process output channels that lack explicit names (i.e., no `emit:` label in the process definition) can be addressed using square bracket index notation. For example, `process.out[0]` accesses the first output channel.

When a channel contains a tuple with multiple elements, you can use additional indices to access specific elements within that tuple, such as `process.out[0][1]` to access the second element of the first channel.

For example, the CAT process from earlier has two output channels with explicit emit labels:

```groovy
output:
tuple val(meta), path("${prefix}"), emit: file_out
path "versions.yml"               , emit: versions
```

You can access these outputs either by name or by index:

- `process.out.file_out` or `process.out[0]` - accesses the first output channel `[meta, file]`
- `process.out.versions` or `process.out[1]` - accesses the second output channel `versions.yml`
- `process.out.file_out[1]` or `process.out[0][1]` - accesses just the file from the first channel tuple

Below you will find examples of a range of different types of assertions that you can apply to module channel outputs.

```groovy
// Process completion status
assert process.success
assert process.exitStatus == 0

// Output channels
assert process.out.my_channel != null
assert process.out.my_channel.size() == 3
assert process.out.my_channel.get(0) == "expected_value"

// For unnamed channels, use index notation
assert process.out[0] != null
assert process.out[0].size() == 3

// Group assertions to see all failures at once
assertAll(
    { assert process.success },
    { assert snapshot(process.out).match() }
)
```

:::note
For more nf-test assertion patterns, see the [nf-test assertions examples documentation](./07_assertions.md).
:::

## Creating a new module with tests

Now we can dive deeper into a step-by-step guide for specifically making an nf-core module.

When creating a new module using `nf-core/tools`, a test file is automatically generated based on the template.

```bash
# Create a new module using nf-core/tools (EXAMPLE)
cd path/to/modules
nf-core modules create seqtk/sample
```

This creates the following module directory structure:

```tree
modules/nf-core/tool/subtool/
├── main.nf              # Process definition
├── meta.yml             # Module metadata
└── tests/               # Testing directory
    ├── main.nf.test     # Test definitions
    └── nextflow.config  # Optional: test-specific config
```

After writing the module under `main.nf` and being ready to write the tests, we can edit the test file (`tests/main.nf.test`) to specify the test's input data via the `input[0]` and `input[1]` channels in the `when` block.

Then in the `test` block, we can define what we want to go into the snapshot, or in other words what to generate to compare subsequent runs against.

In this case we will record the subsampling of an input FASTQ file (`then:`) through generating `md5sums` of all output from the module (where `md5sums` are the default method of recording a file with nf-test).

```groovy
nextflow_process {

    name "Test Process SEQTK_SAMPLE"
    script "../main.nf"
    process "SEQTK_SAMPLE"

    tag "modules"
    tag "modules_nfcore"
    tag "seqtk"
    tag "seqtk/sample"

    test("sarscov2 - fastq") {

        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
                ]
                input[1] = 10  // Number of reads to sample
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

Once we've added the assertions in the `then` block, we run the tests to create a snapshot of the output:

```bash
nf-core modules test seqtk/sample --profile docker
```

This will execute the tests and generate a snapshot file (`tests/main.nf.test.snap`) for validation.

## Stub mode testing

All nf-core modules require at a minimum a `stub:` section to allow 'dry run'-like functionality for pipelines.
This is something that we always want to test for, even if your module also has tests with real data.

Furthermore, in some cases a module may produce output data that is too big for GitHub actions nodes, or the tests take too long or require too much resources.
Therefore, the only option is to test the module in `stub` mode.

To specify a test should be run in `stub` mode, you need to add the `-stub` option as follows:

```main.nf.test
process "MODULE"
config "./nextflow.config"


test("custom evalue") {
  options "-stub"
    when {
        process {
            // test implementation
        }
    }
}
```

:::note
For modules that can never run full tests due to data being too large or requiring too much resources, this option can alternatively be added at the top of `main.nf.test` to have all tests run in stub mode.
:::

## Testing parameter variations

Some modules may require additional parameters added to the test command to successfully run.

These can be specified using a `params` input and an [`ext.args` variable](https://nf-co.re/docs/guidelines/components/modules#configuration-of-extargs-in-tests) within the process scope of the `nextflow.config` file, which exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

If your module requires a `nextflow.config` file to run, create the file in the module's `tests/` directory and add the following code to use parameters defined in the `when` scope of the test.

```bash
touch modules/nf-core/<tool>/<subtool>/tests/nextflow.config
```

```nextflow.config
process {
  withName: 'MODULE' {
    ext.args = params.module_args
  }
}
```

You do not need to modify the contents of this file any further, except for updating the module name to the expected module.

Then import the config into the `main.nf.test` file and supply the params in the `when` section of the test.

```main.nf.test
process "MODULE"
config "./nextflow.config"

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

## Testing chained modules

In some cases, rather than directly linking to pre-made test-data files within the `when:` scope, it may make sense to run an 'upstream module' in your test that produces the required inputs of the module you want to test.

The `setup` method allows you to specify processes or workflows that need to be executed before the primary `when` block.
It serves as a mechanism to prepare the required input data or set up essential steps prior to the primary processing block 'on the fly'.

Within a `setup` block, you can use the `run` method to define and execute multiple dependent processes or workflows.

Here's a basic dummy example of how a setup block looks:

```groovy
nextflow_process {
    name "Test Process MY_MODULE"
    script "../main.nf"
    process "MY_MODULE"

    test("test with setup") {
        setup {
            run("UPSTREAM_MODULE") {
                script "../../upstream/main.nf"
                process {
                    """
                    input[0] = [
                        [ id:'test' ],
                        file('path/to/input.txt')
                    ]
                    """
                }
            }
        }

        when {
            process {
                """
                input[0] = UPSTREAM_MODULE.out.my_output
                """
            }
        }

        then {
            assert process.success
        }
    }
}
```

:::warning
Please keep in mind that changes in processes or workflows executed in the setup method can result in a failed test of a downstream module.
:::

Now let's look at more explicit examples.

### Global `setup` method (for all tests)

A global `setup` method can be defined for all tests within a `nextflow_process` definition.

The `setup` is applied to multiple `test` cases, ensuring a consistent setup for each test.

This approach is useful when multiple tests share the same setup requirements.

```groovy
nextflow_process {

    name "Test Process ABRICATE_SUMMARY"
    script "../main.nf"
    process "ABRICATE_SUMMARY"
    config "./nextflow.config"

    setup {
        run("ABRICATE_RUN") {
            script "../../run/main.nf"
            process {
                """
                input[0] =  Channel.fromList([
                    tuple([ id:'test1', single_end:false ], // meta map
                        file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)),
                    tuple([ id:'test2', single_end:false ],
                        file(params.modules_testdata_base_path + 'genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz', checkIfExists: true))
                ])
                """
            }
        }
    }

    test("first test") {
        when {
            process {
                """
                input[0] = ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
                """
            }
        }
        then {
            assert process.success
            assert snapshot(process.out).match()
        }
    }

    test("second test") {
        when {
            process {
                """
                input[0] = ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
                """
            }
        }
        then {
            assert process.success
            assert snapshot(process.out).match()
        }
    }
}
```

### Local `setup` method (for a single test)

Alternatively a local `setup` method can be defined for a single test within a `nextflow_process` definition.

The `setup` is applied to the specific test, ensuring a consistent setup for that test.

This approach is useful when a test requires a specific setup that is different from the global setup.

Here's an example for `abricate/summary`, which requires output from `abricate/run`:

```groovy
nextflow_process {

    name "Test Process ABRICATE_SUMMARY"
    script "../main.nf"
    process "ABRICATE_SUMMARY"

    tag "modules"
    tag "modules_nfcore"
    tag "abricate"
    tag "abricate/summary"

    test("bacteroides_fragilis - genome_fna_gz") {
        setup {
            run("ABRICATE_RUN") {
                script "../../run/main.nf"
                process {
                    """
                    input[0] = Channel.fromList([
                        tuple([ id:'test1', single_end:false ], // meta map
                            file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)),
                        tuple([ id:'test2', single_end:false ],
                            file(params.modules_testdata_base_path + 'genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz', checkIfExists: true))
                    ])
                    """
                }
            }
        }

        when {
            process {
                """
                // Collect reports from ABRICATE_RUN, create a new meta map, and provide it as input
                input[0] = ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

Run the tests:

```bash
nf-core modules test abricate/summary --profile docker
```

This will execute the test with the chained module setup, running `ABRICATE_RUN` first to generate the required input, then testing `ABRICATE_SUMMARY` with that output.

### Aliasing dependencies

If you need to run the same setup process multiple times for the same test but for different files (such as unzipping multiple different files), you can set an `alias` for the process:

```groovy
nextflow_process {

    // ...

    setup {

        run("UNTAR", alias: "UNTAR1") {
            script "modules/nf-core/untar/main.nf"
            process {
            """
            input[0] = Channel.fromList(...)
            """
            }
        }

        run("UNTAR", alias: "UNTAR2") {
            script "modules/nf-core/untar/main.nf"
            process {
            """
            input[0] = Channel.fromList(...)
            """
            }
        }

        run("UNTAR", alias: "UNTAR3") {
            script "modules/nf-core/untar/main.nf"
            process {
            """
            input[0] = Channel.fromList(...)
            """
            }
        }
    }

    test("Test with three different inputs") {
        when {
            process {
                """
                input[0] = UNTAR1.out.untar.map{ it[1] }
                input[1] = UNTAR2.out.untar.map{ it[1] }
                input[2] = UNTAR3.out.untar.map{ it[1] }
                """
            }
        }

        then {
            // ...
        }
    }
}
```

## Running and maintaining tests

### Updating module snapshots

Whenever a module is updated that results in changes to output (e.g., due to version bumps of the tool itself), you will need to update snapshots:

```bash
nf-core modules test abricate/summary --profile docker --update
```

You will see the following warning at the start of the test run:

```bash
│ 🚀 nf-test 0.9.0                                                                                       │
│ https://www.nf-test.com                                                                                │
│ (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr                                                   │
│                                                                                                        │
│ Load .nf-test/plugins/nft-bam/0.5.0/nft-bam-0.5.0.jar                                                  │
│ Load .nf-test/plugins/nft-compress/0.1.0/nft-compress-0.1.0.jar                                        │
│ Load .nf-test/plugins/nft-vcf/1.0.7/nft-vcf-1.0.7.jar                                                  │
│ Load .nf-test/plugins/nft-csv/0.1.0/nft-csv-0.1.0.jar                                                  │
│ Load .nf-test/plugins/nft-utils/0.0.3/nft-utils-0.0.3.jar                                              │
│ Load .nf-test/plugins/nft-fastq/0.0.1/nft-fastq-0.0.1.jar                                              │
│ Load .nf-test/plugins/nft-anndata/0.1.0/nft-anndata-0.1.0.jar                                          │
│ Warning: every snapshot that fails during this test run is re-record.
```

If you are using nf-core tools for executing the tests, once the test passes, the snapshot will be updated and the test(s) will be re-run for you to verify the snapshot is stable.

If it is not stable, you will need to refine your assertions.

### Running tests with nf-test directly

If you are running with nf-test directly (rather than using `nf-core modules test`), you can use the `--tag` option to specify which module within the repository to test:

```bash
nf-test test --profile docker --tag abricate/summary

# or specify test path
nf-test test --profile docker modules/nf-core/abricate/summary/tests/main.nf.test

# update snapshots
nf-test test --profile docker modules/nf-core/abricate/summary/tests/main.nf.test --update-snapshot
```

This will run the tests for the module and display the results, including any failures or snapshot mismatches.

:::note
The `nf-test test` command runs the test only once compared to the `nf-core modules test` command which runs the test twice to confirm snapshot stability.
:::

## Module testing principles

- **Prefer automated MD5 checksums using Snapshots** for output verification when possible, then file content checks, then existence checks as fallbacks
- **Test both regular process and stub modes** to verify functionality and stub outputs
- **Use appropriate test data** from the **nf-core test-datasets** repository
- **Minimal viable tests** that cover the core functionality without excessive complexity

:::note
**Key points:**

For detailed testing guidelines, see the [nf-core modules testing specifications](https://nf-co.re/docs/guidelines/components/modules#testing).
:::

## Next steps

Continue to [Testing Subworkflows](./04_testing_subworkflows.md) to learn about testing more complex multi-module components.
