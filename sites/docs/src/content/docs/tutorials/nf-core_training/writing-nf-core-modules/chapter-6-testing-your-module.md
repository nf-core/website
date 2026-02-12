---
title: "Chapter 6: Testing your module"
subtitle: "How to write and generate the files needed for testing your module"
shortTitle: "Chapter 6: Testing"
---

## Introduction

In this chapter we will describe how nf-core modules are used for unit testing the modules with the `nf-test` framework.

Once we have the Nextflow files written, we want to make sure the module actually works as intended.

:::warning
Setting up of and debugging of tests is a whole topic in it's own right, and thus could have a dedicated workshop to it.
This chapter provides aims to give you a sufficient introduction to you get started, as tests are _required_ for nf-core modules.
While this step may take a bit of time to get right, don't be discouraged if you don't get it right the first time.
:::

## The `main.nf.test` file

In this file we have 3 main sections

- File paths and tags
- Test block
  - When block
  - Then (assertion) block
- Setup block (optional)

:::info{title="Click here to see full 'raw' file example" collapse}

The boilerplate TODO comments have been removed for readability.

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

### The file paths and tags

This section is initially propagated with the basic information for the test file.

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

It includes an 'umbrella' name for the collection of tests based on the name of the module, the location of the script, and name of the module's process.
It also includes automatically generated tags that are used by nf-test to run the tests in different configurations.

In most cases you do not need to edit this section.

You may however need to extend this section based on how you set up the test later in the file (see below).
These extensions are typically to either add a path to an optional `nextflow.config`, or add additional tags of modules used in the setup block (for both see below).
The latter is to ensure if an 'upstream' module to your module changes, it does not break the test of your new module.

### The `test` block

This section is the main part of the test file.
Here you specify inputs to the test, and how to compare the differences of the module outputs between two test runs.

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

This section is the part of the file you will edit the most.
The test block declaration consists of three main sections:

- The title of the test
- The when block (defining inputs)
- The then block (defining outputs)

The boilerplate template comes with example code for a single input file as input, and tests both that the module successfully passes with no errors, and generates a 'snapshot' of md5sums of the contents of all output channels of the module.
In addition, the boilerplate template also includes a stub runs test.
All nf-core modules require a stub-run test, however you do not need to change this except for the test title, and the inputs (the 'then' block) so they match that the first test.

For writing the test, you need to follow the following steps.

#### `when` block

First, you must update the name of the test to make it distinct.
Typically, at a minimum this will consist of the organism of the test data you will use, and the file format of the primary input file.
You can then also provide an additional keyword or short sentence of what configuration the given test will be testing, for example:

```nextflow
test("sarscov2 - fastq - pairedend")
```

or

```nextflow
test("sarscov2 - fastq - bakta annotation input")
```

Second, within the 'when' block you need to specify each input channel to your module using the notation of `input[index]`, i.e. `input[0]` for the first input channel, `input[1]`, for the second, and so forth.
This 'when' block can be filled with standard Nextflow code - e.g. using a `Channel` factory to create a channel, use operators, and so on, as you normally would in an Nextflow pipeline.
Make sure that the input into each channel `input[]` variable matches the channel restructure of your module (e.g. if there is a meta map)

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

All file inputs should be referred to by utilising files hosted on the [nf-core/test-datasets GitHub repository](https://github.com/nf-core/test-datasets).
You can then load these with the syntax as in the example above (i.e., with `params.modules_testdata_base_path` + the file path as within the nf-core/testdatasets repository), with the `checkIfExists` option.

For optional input channels, these can just be given to the relevant `input` variable with `[]`, e.g. `input[3] = []`, or if it requires a meta map `input[3] = [[],[]]`.
To refer to the output modules run in a setup block (see below), you can refer to these as you would in a Nextflow pipeline, e.g. `PROCESS_NAME.out.foo`.

#### `then` block

Third, once you have completed the input declarations, you can move onto the 'then' block.
This is where you will need to write 'assertions', i.e., telling nf-test what to compare between test runs.
By default, the example in the boilerplate template code will check the contents of every single output file of the channel.
However, very often you may find some tools produce variable output files between runs (a common one is logs which include timestamps).
You can expand the `snapshot()` function to explicitly refer to specific channels, and test each one in different ways.

```nextflow
then {
    assertAll(
        { assert process.success },
        { assert snapshot(
            process.out.metrics,
            process.out.versions,
            file(process.out.qc_report[0][1]).name,
        ).match()
        }
    )
}
```

For example, in the example above you can see that the `metrics` and `versions` output channels do not vary in their md5sum content.
However the `qc_report` does, therefore instead of using the (default) `md5sums` check, we change the 'assertion' so that we compare that name of the files in that channel does not change between runs.

There are many different `assert` methods.
For a more comprehensive list of different nf-test assertions, see the dedicated [nf-core documentation](https://nf-co.re/docs/contributing/nf-test/assertions) page.

We generally recommend the to test files with the following methods in order of preference:

1. 😃 `md5sum`s
2. 🙂 String contents within a file
3. 🫤 File name / file existence

If you have to vary the types of tests per channel, make sure to always include testing of the `process.out.versions` channel!

Finally, once you've completed both the `then` and `when` blocks, you can copy the structure of this first test to each subsequent test case you need to cover - updating the name test and contents.

As a guide, you should try and have as many tests so you test as many configurations as possible, but at a minimum at least enough tests so that you test input files for every input and output channel (mandatory and optional) at least once.

### The `setup` block (optional)

This optional section is where you can specify module(s) to execute _before_ your new module.
The outputs of this upstream module can then be used as input for your new module.
Note that it is not included by default in the boilerplate template code.

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

The location of a setup block can vary.
You can either specify this before all tests, so you can re-use the same output of this upstream module in all subsequent tests.
Alternatively you can place it before the 'when' block of each test within the test block themselves.
In this case the the setup block will only be executed when the given test is executed.

The difference between the setup block and the test block (see above) is that the output of modules in the setup block will _not_ be asserted by the test.

You can fill this block in just the same way as the `test` block, except you must explicitly specify the script path of the upstream module in each setup block.

Otherwise you specify the inputs with the same `input[0]`, `input[1]` etc. channel syntax, and using the URLs to the nf-core test-dataset repository as before.

Note that we generally discourage the use of setup blocks as they increase the runtime of tests.
However they can be useful when a module requires inputs with large file-sizes that are too large for the nf-core/test-datasets repository, or the upstream module is extremely quick.

:::tip{title="Examples" collapse}
Example of a global setup block, the output of which can be reused in every test:

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

Example of setup blocks used for a a specific test:

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

This is an optional file that is _not_ automatically generated with the boilerplate template files.
It sits alongside the `main.nf.test` and `main.nf.test.snap` file within the `tests/` directory.

The purpose of this file is to provide additional settings, e.g. specifying `ext.args` for modules when you want to expand the number of configurations to test.
For example, you could use it if you need to give an optional parameter to produce an additional optional output file.
This can also be use to set up different options for when you need to import a module multiple times (typically in `setup` blocks).

If you need this functionality, you will need to create the file yourself.

```tree {7}
├── environment.yml
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    ├── main.nf.test.snap
    └── nextflow.config
```

You can then specify to use the file within the test by adding this line in your `main.nf.test`:

```nextflow
config './nextflow.config'
```

It should be placed next to the `script` and `tags` declaration. or within each test's scope before the `when` scope (see examples).

:::tip{title="Examples" collapse}
An example of a 'global' nextflow.config file used for all tests:

```nextflow
process {
    ext.args = '--tbi'
}
```

and loaded in the test file

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

or loaded for a specific named test

```nextflow
    test("sarscov2 - paired-end - fastq") {

        config "./nextflow.config"

        setup {
...
```

:::

## The `main.nf.test.snap` file

This file is _not_ automatically generated with the boilerplate template files.
It sits alongside the `main.nf.test` and file within the `tests/` directory, and is automatically generated for you when you run the `main.nf.test` script file with `nf-test` the first time.

```tree {6}
├── environment.yml
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    ├── main.nf.test.snap
    └── nextflow.config
```

You can generate the file by running nf-core tools command:

```bash
nf-core modules test <tool>/<subcommand>
```

and after setting the options via the prompts, the tests will run twice and tell you if the output of each tests were the same or not, in addition to generating the snapshot file (in the first time)

The generated file includes the 'snapshot' output of nf-test.
This snapshot is what will be compared against for each subsequent test run (e.g., comparing the `md5sum` or presence of a string in the output between two given tests).

You never need to manually modify this file, as it will also be automatically updated by nf-test.
If you find variability or an issue in your snapshot, you can re-run with

```bash
nf-core modules test <tool>/<subcommand> --update
```

If you need more information about the test run, you can also provide `--verbose` to see the full nf-test stdout log.

It is important to manually inspect the snapshot file once it has been generated.
You need to make sure you have no unexpected snapshot assertion outputs, even if the tests 'pass' (e.g. empty file `md5sum`, or missing assertions etc.).

If you find issues, and need to inspect the output files of the module itself, we can check the working directory of the module during the test.

To inspect this, in the console when you're running each test, before the name of each test you should get a little 'hash' string:

```bash
│   Test [528b411a] 'candidatus_portiera_aleyrodidarum proteome [fasta]'       │
│ PASSED (7.762s)
```

With this hash string you can change into `.nf-test/tests/<hash string>/work` (make sure to autocomplete with your TAB key to get the full hash) and you can find all the standard Nextflow working directories of the module's test run.
Here, you should change into each directory until you find the one of the module itself, and check the contents of each output file are as expected to ensure your snapshot has been generated correctly.

Note that you may have 'empty' entries in the snapshot file when an optional channel does not emit a file in that given test - however you should double check that is expected for that test.
Furthermore only assertions that are not included in the `snapshot()` function of the `when` block of `main.nf.test` will not be recorded in the snapshot.
Assertions outside the `snapshot()` functions are typically just boolean checks (e.g. existence of a file).
If this type of assertion fails, nf-test will simply fail and be printed to console as an error - not recorded in the snapshot file itself.

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

Once your snapshot is generated and stable, your module is ready for use!

## Summary

In this chapter, we have described how to construct the required and optional files needed for testing your module.

In the next chapter, we will summarise the general development workflow and subsequent contribution to the [nf-core/modules](https://github.com/nf-core/modules) repository.
