---
title: "Testing nf-core/modules"
subtitle: Writing comprehensive tests for nf-core modules with nf-test
weight: 30
---

## Overview

The nf-test framework in the nf-core ecosystem enables comprehensive testing for processes, workflows, and pipelines.

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

- Script paths starting with `./` or `../` are relative to the test script's location.
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

- Each module should contain a `tests/` folder alongside its `main.nf` file
- Test files come with snapshots of component output channels
- Tests verify both functionality and expected outputs
- Support both regular and stub testing modes

## Testing an existing module

Let's examine testing the `bedtools/bamtobed` module, which is a simple standalone module:

```bash
cd path/to/modules
# Run all tests for the module
nf-core modules test bedtools/bamtobed --profile docker
```

This will run all tests for the module and display the results, including any failures or snapshot mismatches.

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

## Creating a new module with tests

When creating a new module using `nf-core/tools`, a test file is automatically generated based on the template.

```bash
# Create a new module using nf-core/tools
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

The generated test file (`tests/main.nf.test`), with positional input channels `input[0]` and `input[1]` provided, will look like this:

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
│ Test Process SEQTK_SAMPLE                                                                              │
│                                                                                                        │
│   Test [a4edc395] 'sarscov2_sample_singleend_fqgz' PASSED (14.223s)                                    │
│   Test [42c1ef08] 'sarscov2_sample_pairedend_fqgz' PASSED (3.353s)                                     │
│   Test [7b327705] 'sarscov2_sample_singlend_fqgz_stub' PASSED (3.042s)                                 │
│   Test [f1c40b60] 'sarscov2_sample_singleend_frac' PASSED (3.326s)                                     │
│                                                                                                        │
│ Test Workflow FASTQ_CONTAM_SEQTK_KRAKEN                                                                │
│                                                                                                        │
│   Test [609e9434] 'sarscov2 - fastq - 25000 - krakendb' PASSED (31.378s)                               │
│   Test [b5fb1cc0] 'sarscov2 - fastq - [12500, 25000] - krakendb' PASSED (7.803s)                       │
│   Test [aac34908] 'sarscov2 - fastq - [12500, 25000] - krakendb -- stub' PASSED (7.262s)               │
│                                                                                                        │
│                                                                                                        │
│ SUCCESS: Executed 7 tests in 70.391s                                                                   │
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
│ Test Process SEQTK_SAMPLE                                                                              │
│                                                                                                        │
│   Test [a4edc395] 'sarscov2_sample_singleend_fqgz' PASSED (3.611s)                                     │
│   Test [42c1ef08] 'sarscov2_sample_pairedend_fqgz' PASSED (3.321s)                                     │
│   Test [7b327705] 'sarscov2_sample_singlend_fqgz_stub' PASSED (3.098s)                                 │
│   Test [f1c40b60] 'sarscov2_sample_singleend_frac' PASSED (3.263s)                                     │
│                                                                                                        │
│ Test Workflow FASTQ_CONTAM_SEQTK_KRAKEN                                                                │
│                                                                                                        │
│   Test [609e9434] 'sarscov2 - fastq - 25000 - krakendb' PASSED (5.934s)                                │
│   Test [b5fb1cc0] 'sarscov2 - fastq - [12500, 25000] - krakendb' PASSED (6.564s)                       │
│   Test [aac34908] 'sarscov2 - fastq - [12500, 25000] - krakendb -- stub' PASSED (7.03s)                │
│                                                                                                        │
│                                                                                                        │
│ SUCCESS: Executed 7 tests in 32.853s                                                                   │
│                                                                                                        │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────╯
INFO     All tests passed!
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

With the config file in place, you can now modify parameters in the `when` block for testing.

```groovy
when {
  params {
    module_args = "--help"
  }
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
```

## Testing chained modules

For modules that depend on the output of other modules, use `setup` blocks to define the dependency. Here's an example for `abricate/summary`, which requires output from `abricate/run`:

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

## Updating module snapshots

When module outputs change (e.g., due to version bumps), you need to update snapshots:

```bash
nf-core modules test abricate/summary --profile docker --update-snapshot
```

:::note
For more nf-test assertion patterns, see the [nf-test assertions examples documentation](./07_assertions.md).
:::

## Next steps

Continue to [Testing Subworkflows](./04_testing_subworkflows.md) to learn about testing more complex multi-module components.
