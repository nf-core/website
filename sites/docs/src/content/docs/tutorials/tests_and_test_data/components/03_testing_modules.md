---
title: "Testing nf-core/modules"
subtitle: Writing comprehensive tests for nf-core modules with nf-test
weight: 30
---

## Overview

The nf-test framework in the nf-core ecosystem enables comprehensive testing for processes, workflows, and pipelines.

This chapter covers the fundamentals of nf-core module testing, from basic syntax to advanced scenarios involving chained modules.

## Module structure

Before diving into testing, let's understand the typical structure of a nf-core module:

```
modules/nf-core/tool/subtool/
├── main.nf              # Process definition
├── meta.yml             # Module metadata
└── tests/               # Testing directory
    ├── main.nf.test     # Test definitions
    └── nextflow.config  # Optional: test-specific config
```

### Basic test syntax

The basic syntax for a process test follows this structure:

```groovy title="main.nf.test"
nextflow_process {
    name "<NAME>"
    script "<PATH/TO/NEXTFLOW_SCRIPT.nf>"
    process "<PROCESS_NAME>"

    test("<TEST_NAME>") {
        when {
        }

        then {
        }
    }
}
```

:::note
**Key points:**

- Script paths starting with `./` or `../` are relative to the test script's location.
  :::

### Essential Assertions

Tests use assertions to verify the expected output of the process specified in the `then` block.

You can specify multiple assertions to be evaluated together in a single test by specifying them within an `assertAll` block.

Nextflow process output channels that lack explicit names (i.e., when no `meta` map is present) can be addressed using square brackets and the corresponding index, for example `process.out[0]` for the first channel and `process.out[0][1]` for the second element of the first channel.

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

## Module testing principles

- **Prefer automated MD5 checksums using Snapshots** for output verification when possible, then file content checks, then existence checks as fallbacks
- **Test both regular process and stub modes** to verify functionality and stub outputs
- **Use appropriate test data** from the **nf-core test-datasets** repository
- **Minimal viable tests** that cover the core functionality without excessive complexity

:::note
**Key points:**

- For detailed testing guidelines, see the [nf-core modules testing specifications](https://nf-co.re/docs/guidelines/components/modules#testing).
  :::

## Creating a new module with tests

When creating a new module using `nf-core/tools`, a test file is automatically generated based on the template.

```bash
# Create a new module using nf-core/tools (EXAMPLE)
cd path/to/modules
nf-core modules create seqtk/sample
```

This creates the following module directory structure:

```
modules/nf-core/seqtk/sample/
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
```

The generated test file (`tests/main.nf.test`) will look like this, once the `input[0]` and `input[1]` channels are defined in the `when` blocks:

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

    test("sarscov2 - fastq - stub") {

        options "-stub"

        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
                ]
                input[1] = 10
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

Run the tests to create a snapshot of the output:

```bash
nf-core modules test seqtk/sample --profile docker
```

This will execute the tests and generate snapshots for validation.

## Testing an existing module

Let's examine testing the `bedtools/bamtobed` module, which is a simple standalone module:

```bash
cd path/to/modules
# Run all tests for the module
nf-core modules test bedtools/bamtobed --profile docker
```

This will run all tests for the module specified in the `main.nf.test` and display the results, including any failures or snapshot mismatches.

```bash


                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.3.2 - https://nf-co.re


INFO     Generating nf-test snapshot
╭──────────────────────────────────────────── nf-test output ────────────────────────────────────────────╮
│                                                                                                        │
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
│                                                                                                        │
│ Test Process BEDTOOLS_BAMTOBED                                                                         │
│                                                                                                        │
│   Test [824188d1] 'sarscov2 - bam' PASSED (3.335s)                                                     │
│   Test [f4f6429b] 'stub' PASSED (3.154s)                                                               │
│                                                                                                        │
│                                                                                                        │
│ SUCCESS: Executed 2 tests in 6.498s                                                                    │
│                                                                                                        │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────╯
INFO     Generating nf-test snapshot again to check stability
╭──────────────────────────────────────────── nf-test output ────────────────────────────────────────────╮
│                                                                                                        │
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
│                                                                                                        │
│ Test Process BEDTOOLS_BAMTOBED                                                                         │
│                                                                                                        │
│   Test [824188d1] 'sarscov2 - bam' PASSED (3.277s)                                                     │
│   Test [f4f6429b] 'stub' PASSED (3.161s)                                                               │
│                                                                                                        │
│                                                                                                        │
│ SUCCESS: Executed 2 tests in 6.446s                                                                    │
│                                                                                                        │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────╯
INFO     All tests passed!
```

### Testing parameter variations

Some modules MAY require additional parameters added to the test command to successfully run.

These can be specified using a params input and an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

If your module requires a `nextflow.config` file to run, create the file to the module’s `tests/` directory and add the following code to use parameters defined in the `when` scope of the test.

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

You do not need to modify the contents of this file any further.

Then add the config to the main.nf.test file and supply the params in the when section of the test.

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

### Choosing Parameter Configuration Methods

**Use `nextflow.config` when:**

- Multiple tests in the same module need the same parameter structure
- You need complex process-specific configurations (memory, cpus, etc.)
- Parameters require process selectors or conditional logic
- You want to maintain consistent configuration across all tests

**Use inline `params` blocks when:**

- Testing different parameter values in individual tests
- You need test-specific parameter overrides
- Parameters are simple and don't require process configuration
- You want to keep test parameters close to the test logic

**Example comparison:**

```groovy
// nextflow.config approach - good for consistent configuration
process {
  withName: 'BLAST_BLASTN' {
    ext.args = params.blast_args
    memory = params.blast_memory ?: '8.GB'
  }
}

// Inline params approach - good for test-specific values
test("custom evalue") {
    when {
        params {
            blast_args = '-evalue 0.001 -max_target_seqs 10'
        }
        process {
            // test implementation
        }
    }
}
```

### Stub mode

When your test data is too big, the tests take too long or require too much resources, you can opt to run your tests in stub mode by adding the following option:

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
This can be added at the top of `main.nf.test` to have all tests run in stub mode or this can also be added to a single test
:::

## Testing chained modules

In some cases, rather than directly linking to pre-made test-data files, it may make sense to run an 'upstream module' in your test to output the required inputs of the module you want to test.

The `setup` method allows you to specify processes or workflows that need to be executed before the primary `when` block.

It serves as a mechanism to prepare the required input data or set up essential steps prior to the primary processing block.

Within the `setup` block, you can use the `run` method to define and execute multiple dependent processes or workflows.

Here's a basic example of how a setup block looks:

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
Please keep in mind that changes in processes or workflows, which are executed in the setup method, can result in a failed test run.
:::

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

A local `setup` method can be defined for a single test within a `nextflow_process` definition.

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

```bash


                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.3.2 - https://nf-co.re


INFO     Generating nf-test snapshot
╭──────────────────────────────────────────── nf-test output ────────────────────────────────────────────╮
│                                                                                                        │
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
│                                                                                                        │
│ Test Process ABRICATE_SUMMARY                                                                          │
│                                                                                                        │
│   Test [fc133477] 'Should run without failures' PASSED (102.456s)                                      │
│                                                                                                        │
│                                                                                                        │
│ SUCCESS: Executed 1 tests in 102.459s                                                                  │
│                                                                                                        │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────╯
INFO     Generating nf-test snapshot again to check stability
╭──────────────────────────────────────────── nf-test output ────────────────────────────────────────────╮
│                                                                                                        │
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
│                                                                                                        │
│ Test Process ABRICATE_SUMMARY                                                                          │
│                                                                                                        │
│   Test [fc133477] 'Should run without failures' PASSED (11.658s)                                       │
│                                                                                                        │
│                                                                                                        │
│ SUCCESS: Executed 1 tests in 11.667s                                                                   │
│                                                                                                        │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────╯
INFO     All tests passed!
```

### Aliasing dependencies

If you need to run the same setup process multiple times for the same test but for different files, you can set an `alias` for the process:

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

## Updating module snapshots

Whenever a module is updated resulting in changes to output (e.g., due to version bumps of the tool itself), you will need to update snapshots:

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

## Testing with nf-test

If you are running with nf-test directly, you can use the `--tag` option to specify which module within the repository to test:

```bash
nf-test test --profile docker --tag abricate/summary

# or specify test path
nf-test test --profile docker modules/nf-core/abricate/summary/tests/main.nf.test

# update snapshots
nf-test test --profile docker modules/nf-core/abricate/summary/tests/main.nf.test --update-snapshot
```

This will run the tests for the module and display the results, including any failures or snapshot mismatches.

```bash
🚀 nf-test 0.9.0
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Load .nf-test/plugins/nft-bam/0.5.0/nft-bam-0.5.0.jar
Load .nf-test/plugins/nft-compress/0.1.0/nft-compress-0.1.0.jar
Load .nf-test/plugins/nft-vcf/1.0.7/nft-vcf-1.0.7.jar
Load .nf-test/plugins/nft-csv/0.1.0/nft-csv-0.1.0.jar
Load .nf-test/plugins/nft-utils/0.0.3/nft-utils-0.0.3.jar
Load .nf-test/plugins/nft-fastq/0.0.1/nft-fastq-0.0.1.jar
Load .nf-test/plugins/nft-anndata/0.1.0/nft-anndata-0.1.0.jar

Test Process ABRICATE_SUMMARY

  Test [fc133477] 'Should run without failures' PASSED (15.286s)


SUCCESS: Executed 1 tests in 15.294s
```

:::note
The `nf-test test` command runs the test only once compared to the `nf-core modules test` command which runs the test twice to confirm snapshot stability.
:::

:::note
For more nf-test assertion patterns, see the [nf-test assertions examples documentation](./07_assertions.md).
:::

## Next steps

Continue to [Testing Subworkflows](./04_testing_subworkflows.md) to learn about testing more complex multi-module components.
