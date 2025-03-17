---
title: "Chapter 6: Testing your module"
subtitle: "How to write and generate the files needed for testing your module"
shortTitle: "Chapter 6: Testing"
---

In this chapter we will describe how nf-core modules used for unit testing the modules with the `nf-test` framework.

## The `main.nf.test` file

In this file we have X main sections

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

```
nextflow_process {

    name "Test Process DREP_COMPARE"
    script "../main.nf"
    process "DREP_COMPARE"

    tag "modules"
    tag "modules_nfcore"
    tag "drep"
    tag "drep/compare"
```

It includes the umbrella name of the tests based on the names of the module, and the location and name of the module script itself.
It also includes automatically generated tags that are used by nf-test to run the tests in different configurations.

In most cases you do not need to edit this section.

You may however have to extend this section based on how you set up the test later in the file (see below).
These extensions typically are to either add a path to a `nextflow.config`, or add additional tags of modules used in the setup block (see below).
The latter is to ensure if an 'upstream' module to your module changes, it does not break the test of your new module.

### The `test` block

This section is the main part of the test file, where you specify inputs to the test, and how to compare the differences of the module outputs between two test runs.

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
The test block declaration consists of three main sections: the title of the test, the when block (defining inputs), and the then block (defining outputs).

The boilerplate template comes with example code for a single input file as input, and tests both that the module successfully passes with no errors, and generate a 'snapshot' of md5sums of the contents of all output channels of the module.
In addition, the boilerplate template also includes a stub runs test.
All nf-core modules require a stub-run test, however you do not need to change this except for the test title, and the inputs (the 'then' block) so they match that the first test.

First, you must update the name of the test to make it distinct.
Typically, at a minimum this will consist of the organism of the test data you will use, and the file format of the primary input file.
You can then also provide an additional keyword or short sentence of what configuration the given test will be testing, for example `test("sarscov2 - fastq - pairedend")`or `test("sarscov2 - fastq - bakta annotation input")`

Within the 'when' block, you then need to specify each input channel to your module using the notation of `input[0]` for the first input channel, `input[1]`, for the second, and so forth.
This 'when' block can be filled with standard Nextflow code - e.g. using a `Channel` factory to create a channel, use operators, and so on, as you normally would in an Nextflow pipeline.
Make sure that the input into each channel `input[]` variable matches the channel restructure of your module (e.g. if there is a metamap)
All file inputs should be referred to by utilising files hosted on the [nf-core/test-datasets GitHub repository](https://github.com/nf-core/test-datasets).
You can then load these with the syntax as in the example above (i.e., with `params.modules_testdata_base_path` + the file path as within the nf-core/testdatasets repository), with the `checkIfExists` option.
For optional input channels, these can just be given to the relevant `input` variable with `[]`, e.g. `input[3] = []`, or if it requires a meta map `input[3] = [[:],[]]`.
To refer to the output modules run in a setup block (see below), you can refer to these as you would in a Nextflow pipeline, e.g. `PROCESS_NAME.out.foo`.

Once you have completed the input declarations, you can move onto the 'then' block.
This is where you will need to write 'assertions', i.e., telling nf-test what to compare for.
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
However the `qc_report` does, therefore instead of `md5sums` we change the 'assertion' so that the name of the files in that channel does not change between runs.
For a more comprehensive list of different nf-test assertions, see the dedicated [nf-core documenation](https://nf-co.re/docs/contributing/nf-test/assertions) page.

If you have to vary the types of tests per channel, make sure to always include testing of the `process.out.versions` channel!

Once you've completed both the `then` and `when` blocks, you can copy the structure and update the name test and contents for every subsequent tests.
As a guide, you should try and have as many tests so you test as many configurations as possible, but at a minimum at least a test so you test input files for all input and output channels (mandatory and optional) at least once.

### The `setup` block (optional)

This optional section is where you can specify module(s) to execute _before_ your new module itself, to use the outputs of this upstream module as input for your module.
It is not included by default in the boilerplate template code.

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

The difference between this and the test block (see below) is that output of modules in the setup block will not be evaluated by the test.

You can fill this block in just the same way as the `test` block, except you explicitly specify the script path of the upstream module.

Otherwise you specify the inputs in the same `input[0]`, `input[1]` etc channel syntax, and using the URLs to the nf-core test-dataset repository.

Note that we generally discourage use of setup blocks unless necessary to reduce runtime of tests.
However it can be useful when a module requires inputs with large file-sizes, or the upstream module is extremely quick.

## The `nextflow.config` file (optional)

This is an optional file that is _not_ automatically generated with the boilerplate template files.
It sits alongside the `main.nf.test` and `main.nf.test.snap` file within the `tests/` directory.

The purpose of this file is to provide additional settings, e.g. specifying `ext.args` for modules when you want to expand the number of tests are run.
This could be for example if you need to give an optional parameter to produce an additional optional output file.
This can also be use to set up different options for when you need to import a module multiple times (typically in `setup` blocks).

If you create such a file, this file you will need to specify the test to use this file, by adding this line:

```nextflow
config './nextflow.config'
```

Either next to the `script` and `tags` declaration. or within each test's scope before the `when` scope.

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

You can generate the file by running nf-core tools command:

```bash
nf-core modules test <tool>/<subcommand>
```

and after setting the options via the dialogue, the tests will run twice and tell you if the results of the tests were the same or not, in addition to generating the snapshot file (in the first time)

The generated file includes the 'snapshot' output results of nf-test, that are compared against against each subsequent test (e.g., compared the `md5sum` or presence of a string in the output between two tests runs).

You never need to manually modify this file, as it will also be automatically updated by nf-test.

However it is important to manually inspect this file once it has been generated that you have no unexpected snapshot assertion outputs (e.g. empty file `md5sum`, or missing assertions etc.

To inspect this, in the console when you're running each test, before the name of each test you should get a little 'hash' string:

```bash
│   Test [528b411a] 'candidatus_portiera_aleyrodidarum proteome [fasta]'       │
│ PASSED (7.762s)
```

With this hash string you can change into `.nf-test/tests/528b411a<rest of the string>/work` and you can find a standard Nextflow working directory of the module's test run.
Within here you should change into each directory until you find the one of the module itself, and check the contents of each output file are as expected to ensure your snapshot has been generated correctly.

Note that you may have 'empty' entries when an optional channel does not emit a file in that given test - however you should double check that is expected for that test.
Furthermore only assertions that are not included in the snapshot function will not be recorded in the snapshot (typically just evaluated by a boolean, if failed will be printed to console as an error).

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

Once your snapshot is generated, your module is ready for use!

## Summary

In this chapter, we have described how to construct the required and optional files needed for testing your module.

In the next chapter, we will summarise the general development workflow and subsequent contribution to the [nf-core/modules](https://github.com/nf-core/modules) repository.
