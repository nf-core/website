---
title: "Chapter 6: Testing your module"
subtitle: "How to test an nf-core module"
shortTitle: "Chapter 6: Testing modules"
---

This chapter describes how to unit test an nf-core module with the [nf-test](https://www.nf-test.com/) framework.

:::warning
Setting up and debugging tests is a topic of its own.
This chapter covers the basics to get you started.
Tests are required for all nf-core modules, so do not be discouraged if you need multiple attempts to get them right.
:::

## The `main.nf.test` file

The `main.nf.test` file has three main sections:

- File paths and tags.
- Test block, made up of a `when` block (inputs) and a `then` block (assertions).
- Setup block (optional).

:::info{title="Click here to see full 'raw' file example" collapse}

The boilerplate `TODO` comments have been removed for readability.

```nextflow
nextflow_process {

    name "Test Process DREP_COMPARE"
    script "../main.nf"
    process "DREP_COMPARE"

    tag "modules"
    tag "modules_nfcore"
    tag "drep"
    tag "drep/compare"

    test("sarscov2 - bam") {

        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true),
                ]
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

    test("sarscov2 - bam - stub") {

        options "-stub"

        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true),
                ]
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

:::

### File paths and tags

This section is pre-filled with the module's metadata.

```nextflow
nextflow_process {

    name "Test Process DREP_COMPARE"
    script "../main.nf"
    process "DREP_COMPARE"

    tag "modules"
    tag "modules_nfcore"
    tag "drep"
    tag "drep/compare"
```

It contains:

- An umbrella test name derived from the module name.
- The script location and process name.
- Auto-generated tags that nf-test uses to group tests.

You rarely need to edit this section. You may need to extend it to:

- Add a path to an optional `nextflow.config` (see below).
- Add tags for modules used in a `setup` block, so that upstream module changes do not silently break your test.

### The `test` block

The test block declares inputs and assertions.
You will spend most of your time here.

```nextflow
test("sarscov2 - bam") {
    when {
        process {
            """
            input[0] = [
                [ id:'test', single_end:false ], // meta map
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true),
            ]
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
```

A test block has three parts:

- The test title.
- The `when` block (inputs).
- The `then` block (assertions).

The boilerplate provides one test with a single input file plus a stub-run test.
All nf-core modules need a stub-run test.
You only need to update its title and inputs to match the first test.

#### `when` block

Update the test title to make it distinct.
At a minimum include the organism of the test data and the file format of the primary input, for example:

```nextflow
test("sarscov2 - fastq - pairedend")
```

or

```nextflow
test("sarscov2 - fastq - bakta annotation input")
```

In the `when` block, declare each input channel using `input[index]` notation (`input[0]` for the first, `input[1]` for the second, and so on).
The block accepts standard Nextflow code — channel factories, operators, and so on — as you would use in a pipeline.
Make sure each `input[]` matches the channel structure of your module, including any meta map:

```nextflow
when {
    process {
        """
        ch_samplesheet = Channel.of([
            [ id:'test' ],
            file(params.modules_testdata_base_path + 'genomics/homo_sapiens/array_expression/GSE38751.csv', checkIfExists: true)
            ]
        )
        input[0] = ch_samplesheet.join(UNTAR.out.untar)
        input[1] = [[],[]]
        """
    }
}
```

Key rules:

- Reference all input files via the [nf-core/test-datasets repository](https://github.com/nf-core/test-datasets), using `params.modules_testdata_base_path` with `checkIfExists`.
- For optional input channels, pass `[]`, or `[[],[]]` if a meta map is required.
- To reference output from a `setup` block, use `PROCESS_NAME.out.foo`.

#### `then` block

The `then` block holds assertions — what nf-test compares between runs. The boilerplate checks the contents of every output channel by default.

Some tools produce variable outputs (logs with timestamps, for example). You can expand `snapshot()` to reference specific channels and assert them in different ways:

```nextflow
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            process.out.metrics,
           process.out.findAll { key, val -> key.startsWith('versions') }
            file(process.out.qc_report[0][1]).name,
        ).match()
        }
    )
}
```

Here, `metrics` have stable md5 sums and use the default snapshot check.
The `qc_report` varies, so we only assert that the filename is unchanged.
We use a special assertion to print the versions so they are readable in the test output file.

For a full list of nf-test assertions, see the [nf-core documentation on assertions](https://nf-co.re/docs/developing/testing/assertions).

Test outputs using these methods, in order of preference:

1. `md5sums`.
2. String contents within a file.
3. File name or existence.

If you vary assertion types per channel, always include `process.out.versions`.

Once you have completed one test, copy its structure for each additional test case. Write as many tests as needed to cover every configuration, and at a minimum test every mandatory and optional input and output channel once.

### The `setup` block (optional)

The `setup` block runs one or more modules before your new module. The outputs of these upstream modules can then feed into your test. The setup block is not included by default.

```nextflow
    setup {
        run ("UNTAR") {
            script "../../../untar/main.nf"
            process {
                """
                input[0] = [[],file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/db/kraken2_bracken.tar.gz', checkIfExists: true)]
                """
            }
        }
    }
```

Placement options:

- Before all tests, to reuse the upstream output across every test.
- Inside a test block, to run the setup only for that test.

Outputs produced in a setup block are not asserted against.

Fill in the setup block like a test block, but include the script path of the upstream module. Use the same `input[0]`, `input[1]` syntax.

We discourage setup blocks because they increase runtime. They are useful when the module needs inputs too large for the nf-core/test-datasets repository, or when the upstream module runs quickly.

:::tip{title="Examples" collapse}
A global setup block, reused in every test:

```nextflow
nextflow_process {

    name "Test Process ADAPTERREMOVALFIXPREFIX"
    script "../main.nf"
    process "ADAPTERREMOVALFIXPREFIX"

    tag "modules"
    tag "modules_nfcore"
    tag "adapterremoval"
    tag "adapterremovalfixprefix"

    setup {
        run("ADAPTERREMOVAL") {
            script "../../adapterremoval/main.nf"
            config "./nextflow.config"
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ], // meta map
                    [
                        file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
                        file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true)
                    ]
                ]
                input[1] = []
                """
            }
        }
    }

    test("paired-end - sarscov2 - [fastq]") {

        when {

            process {
                """
                input[0] = ADAPTERREMOVAL.out.collapsed
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
...
```

A setup block scoped to a specific test:

```nextflow
nextflow_process {

    name "Test Process BCFTOOLS_PLUGINIMPUTEINFO"
    script "../main.nf"
    process "BCFTOOLS_PLUGINIMPUTEINFO"

    tag "modules"
    tag "modules_nfcore"
    tag "bcftools"
    tag "bcftools/pluginimputeinfo"
    tag "bcftools/plugintag2tag"

    test("sarscov2 - [vcf, tbi], [], []") {

        config "./nextflow.config"

        setup {
            run("BCFTOOLS_PLUGINTAG2TAG") {
                script "../../plugintag2tag/main.nf"
                process {
                    """
                    input[0] = [
                        [ id:'out', single_end:false ], // meta map
                        file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/vcf/test.vcf.gz', checkIfExists: true),
                        file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/vcf/test.vcf.gz.tbi', checkIfExists: true)
                    ]
                    input[1] = []
                    input[2] = []
                    """
                }
            }
        }

        when {
            process {
                """
                input[0] = BCFTOOLS_PLUGINTAG2TAG.out.vcf.join(BCFTOOLS_PLUGINTAG2TAG.out.tbi)
                input[1] = []
                input[2] = []
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert process.success },
                { assert snapshot(
                    process.out.vcf,
                    process.out.versions
                ).match() }
            )
        }

    }
```

:::

## The `nextflow.config` file (optional)

The optional `nextflow.config` file is not generated by the boilerplate. It sits alongside `main.nf.test` and `main.nf.test.snap` in the `tests/` directory.

Use it to provide extra settings, for example `ext.args` for test configurations, optional parameters that trigger extra outputs, or per-import settings when a module is imported more than once in setup blocks.

Create the file yourself when you need it:

```tree {7}
├── environment.yml
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    ├── main.nf.test.snap
    └── nextflow.config
```

Reference it in the test file:

```nextflow
config './nextflow.config'
```

Place the line next to the `script` and `tags` declarations, or inside a test block before its `when` scope.

:::tip{title="Examples" collapse}
A global `nextflow.config` for all tests:

```nextflow
process {
    ext.args = '--tbi'
}
```

Loaded at the top of the test file:

```nextflow
nextflow_process {

    name "Test Process BCFTOOLS_INDEX"
    script "../main.nf"
    config "./nextflow.config"
    process "BCFTOOLS_INDEX"

    tag "modules"
    tag "modules_nfcore"
    tag "bcftools"
    tag "bcftools/index"

...

```

Or loaded inside a single test:

```nextflow
    test("sarscov2 - paired-end - fastq") {

        config "./nextflow.config"

        setup {
...
```

:::

## The `main.nf.test.snap` file

The `main.nf.test.snap` file is not part of the boilerplate. nf-test generates and updates it automatically the first time you run the test.

```tree {6}
├── environment.yml
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    ├── main.nf.test.snap
    └── nextflow.config
```

Generate it with:

```bash
nf-core modules test <tool>/<subcommand>
```

Set the prompts, and the tests run twice. nf-test reports whether the outputs match and creates the snapshot on the first run.

The snapshot records what every later test run is compared against, such as md5 sums or string matches.

You never need to edit this file manually. If you need to regenerate it:

```bash
nf-core modules test <tool>/<subcommand> --update
```

Add `--verbose` to see the full nf-test output.

Inspect the snapshot after generation. Even if tests pass, look for unexpected values such as empty file md5 sums or missing assertions.

To inspect a module's working directory, note the hash shown next to each test name:

```bash
│   Test [528b411a] 'candidatus_portiera_aleyrodidarum proteome [fasta]'       │
│ PASSED (7.762s)
```

Change into `.nf-test/tests/<hash>/work` (use TAB to autocomplete) to find the Nextflow working directories for that test. Each subdirectory contains the outputs of one process run, so you can verify the contents of every output file.

Empty entries in the snapshot may be legitimate when an optional channel does not emit. Always verify this is expected for the test.

Assertions outside the `snapshot()` function are not recorded in the snapshot. These are usually boolean checks such as file existence, and nf-test prints failures to the console instead.

:::tip{title="Examples" collapse}

```nextflow
{
    "sarscov2 nanopore [fastq_gz]": {
        "content": [
            [
                [
                    {
                        "id": "test"
                    },
                    [
                        "test.fastq_read-assignment-distributions.tsv:md5,c5c52bef375dee47407b3711b147b61d",
                        "test.fastq_rel-abundance.tsv:md5,6fb86cb4103ae5064bb2b7b43e9c9c24"
                    ]
                ]
            ],
            [
                [
                    {
                        "id": "test"
                    },
                    "test.fastq_read-assignment-distributions.tsv:md5,c5c52bef375dee47407b3711b147b61d"
                ]
            ],
            [
                [
                    {
                        "id": "test"
                    },
                    "test.fastq_emu_alignments.sam:md5,09f832656b47f44ff22f65ff7ac87824"
                ]
            ],
            [
                "versions.yml:md5,5227d63df8748de6a8b0c3f212a79512"
            ]
        ],
        "meta": {
            "nf-test": "0.9.2",
            "nextflow": "24.10.4"
        },
        "timestamp": "2025-02-17T15:48:16.876032268"
    },
    "sarscov2 nanopore [fastq_gz] - stub": {
        "content": [
            {
                "0": [
                    [
                        {
                            "id": "test"
                        },
                        "test.tsv:md5,d41d8cd98f00b204e9800998ecf8427e"
                    ]
                ],
                "1": [

                ],
                "2": [
                    [
                        {
                            "id": "test"
                        },
                        "test.sam:md5,d41d8cd98f00b204e9800998ecf8427e"
                    ]
                ],
                "3": [
                    [
                        {
                            "id": "test"
                        },
                        "test.fasta:md5,d41d8cd98f00b204e9800998ecf8427e"
                    ]
                ],
                "4": [
                    "versions.yml:md5,5227d63df8748de6a8b0c3f212a79512"
                ],
                "assignment_report": [

                ],
                "report": [
                    [
                        {
                            "id": "test"
                        },
                        "test.tsv:md5,d41d8cd98f00b204e9800998ecf8427e"
                    ]
                ],
                "samfile": [
                    [
                        {
                            "id": "test"
                        },
                        "test.sam:md5,d41d8cd98f00b204e9800998ecf8427e"
                    ]
                ],
                "unclassified_fa": [
                    [
                        {
                            "id": "test"
                        },
                        "test.fasta:md5,d41d8cd98f00b204e9800998ecf8427e"
                    ]
                ],
                "versions": [
                    "versions.yml:md5,5227d63df8748de6a8b0c3f212a79512"
                ]
            }
        ],
        "meta": {
            "nf-test": "0.9.2",
            "nextflow": "24.10.4"
        },
        "timestamp": "2025-02-17T15:48:24.167181483"
    }
}
```

:::

Once your snapshot is stable, your module is ready for use.

The next chapter summarises the development workflow and the steps for contributing to the [nf-core/modules](https://github.com/nf-core/modules) repository.
