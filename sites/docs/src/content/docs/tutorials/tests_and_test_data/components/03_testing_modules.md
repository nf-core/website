---
title: "Testing Modules"
subtitle: Writing comprehensive tests for nf-core modules with nf-test
weight: 30
---

## Overview

The nf-test framework in nf-core ecosystem enables comprehensive testing process/workflow/pipeline.

This chapter covers the fundamentals of nf-core module testing, from basic syntax to advanced scenarios involving chained modules.

### Basic test syntax

The basic syntax for a process test follows this structure:

```groovy
nextflow_process {
    name "<NAME>"
    script "<PATH/TO/NEXTFLOW_SCRIPT.nf>"
    process "<PROCESS_NAME>"

    test("<TEST_NAME>") {
        // Test implementation
    }
}
```

:::note
**Key points:**

- Script paths starting with `./` or `../` are relative to the test script location
- Use relative paths to reference files within the same directory or parent directories
  :::

### Essential Assertions

Process tests commonly use these assertions:

```groovy
// Process status
assert process.success
assert process.exitStatus == 0

// Output channels
assert process.out.my_channel != null
assert process.out.my_channel.size() == 3
assert process.out.my_channel.get(0) == "expected_value"

// For unnamed channels, use index notation
assert process.out[0] != null
assert process.out[0].size() == 3
```

## Module testing principles

Following the [nf-core testing guidance](https://nf-co.re/docs/tutorials/tests_and_test_data/nf-test_writing_tests), each nf-core module should include comprehensive tests that:

- Each module should contain a `tests/` folder alongside its `main.nf` file
- Test files come with snapshots of component output channels
- Tests verify both functionality and expected outputs
- Support both regular and stub testing modes

## Creating a new module with tests

When creating a new module using nf-core tools, a test file is automatically generated based on the template.

```bash
# Create a new module using nf-core tools
cd path/to/modules
nf-core modules create seqtk/sample
```

This creates the following module structure:

```
modules/nf-core/seqtk/sample/
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    └── tags.yml
```

The generated test file (`tests/main.nf.test`), with inputs[0/1] provided will look like this:

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

After providing the appropriate test data, run the tests to create a snapshot of the output:

```bash
nf-core modules test seqtk/sample --profile docker
```

### Adding parameters to tests

For modules requiring additional parameters, you can create a `<filename>.config` file in the `tests/` directory:

```bash
# Create config file for parameter testing
touch modules/nf-core/seqtk/sample/tests/nextflow.config
```

Add the configuration:

```groovy
process {
    withName: 'SEQTK_SAMPLE' {
        ext.args = params.module_args
    }
}
```

Then apply the config in your test as shown below:

```groovy
nextflow_process {
    name "Test Process SEQTK_SAMPLE"
    script "../main.nf"
    process "SEQTK_SAMPLE"
    config "./nextflow.config"
}
```

With the config file in place, you can now run the tests with modifying the parameter in the `when` block.

```groovy
when {
  params {
    module_args = "--help"
  }
  process {
    """
    input[0] = [
      [ id:'test', single_end:false ], // meta map
      file(params.modules_testdata_base_path + 'genomics/fastq/test_1.fastq.gz')
    ]
    """
  }
}
```

## Testing an existing module

Let's examine testing the `bedtools/bamtobed` module, which is a simple standalone module:

```bash
# Run all tests for the module
nf-core modules test bedtools/bamtobed --profile docker
```

```bash
OUTPUT
```

## Testing chained modules

For modules that depend on outputs from other modules, use the [setup](https://www.nf-test.com/docs/testcases/setup/) method. Here's an example for `abricate/summary`, which requires output from `abricate/run`:

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

```bash
OUTPUT
```

## Updating module snapshots

When module outputs change (e.g., due to version bumps), you need to update snapshots:

```bash
nf-core modules test abricate/summary --profile docker --update-snapshot
```

:::note
For more nf-test assertion patterns, see the [nf-test assertions examples documentation](07_assertions.md).
:::

## Next steps

Continue to [Testing Subworkflows](./04_testing_subworkflows.md) to learn about testing more complex multi-module components.
