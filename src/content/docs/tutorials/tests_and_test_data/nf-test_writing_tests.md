---
title: 'nf-test: Writing tests'
subtitle: Guidelines for writing nf-test tests
shortTitle: Writing nf-test tests
---

## Philosophy of nf-tests

- Each component contains a `tests/` folder beside the `main.nf` of the component itself, containing the test files.
- Test files come with a [snapshot](https://code.askimed.com/nf-test/docs/assertions/snapshots/) of component output channels.

## nf-test guidelines for a simple un-chained module

- Some modules MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

If your module requires a a `nextflow.config` file to run, create the file to the module's `tests/` directory and add the additional parameters there.

```bash
touch modules/nf-core/<tool>/<subtool>/tests/nextflow.config
```

Then add the path to the `main.nf.test` file.

```groovy title="main.nf.test"
process "MODULE"
config "./nextflow.config"
```

- When your test data is too big, the tests take too long or require too much resources, you can opt to run your tests in stub mode by adding the following option:

```groovy title="main.nf.test"
options "-stub"
```

:::note
this can be added at the top of `main.nf.test` to have all tests run in stub mode or this can also be added to a single test
:::

- You can find examples of different nf-tests assertions on [this tutorial](/docs/contributing/nf-test/assertions).

## nf-test guidelines for a chained module

- For modules that involve running more than one process to generate required test-data (aka chained modules), nf-test provides a [setup](https://code.askimed.com/nf-test/docs/testcases/setup/) method.

- For example, the module `abricate/summary` requires the process `abricate/run` to be run prior and takes its output as input. The `setup` method is to be declared before the primary `when` block in the test file as shown below:

```groovy title="main.nf.test"
setup {

            run("ABRICATE_RUN") {
                script "../../run/main.nf"
                process {
                    """
                    input[0] =  Channel.fromList([
                        tuple([ id:'test1', single_end:false ], // meta map
                            file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)),
                        tuple([ id:'test2', single_end:false ],
                            file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true))
                    ])
                    """
                }
            }
        }
```

:::note
The setup method can run more than one process each enclosed in their own `run` block
:::

- Then, the output of setup process/es can be provided as input in the `process` section of `when` block

```groovy title="main.nf.test"
input[0] = ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
```

- Next, in the `then` block we can write our assertions that are used to verify the test. A test can have multiple assertions but, we recommend enclosing all assertions in a `assertAll()` block as shown below:

```groovy title="main.nf.test"
assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
          )
```

- the `main.nf.test` file for chained modules will finally look as shown below:

```groovy title="main.nf.test"
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
                                    file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)),
                                tuple([ id:'test2', single_end:false ],
                                    file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true))
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

## nf-test guidelines for pipelines

Pipeline level tests can facilitate more reliable and reproducible pipelines by ensuring the pipeline produces identical results with every run. You must add pipeline tests that work with `-profile test` and you should reuse this profile within one nf-test.

### Pipeline nf-test overview and structure

Within the base directory of the repository, there is a configuration file for nf-test, named `nf-test.config`. This will set the following options:

- Set the testsDir to the base of the repository so it includes all files
- Set the default profile(s) for nf-test to include `test` (this can be overridden on the command line)
- Add an additional configuration file specific for nf-test located in `tests/nextflow.config`

Within the nf-test specific configuration file, you can add specific requirements for running nf-tests but this should not include parameters or options as these should be available in all contexts.

Pipeline level tests are located in the `tests` directory. Each nf-test file must contain a single pipeline test which tests for a single scenario. Although nf-test supports multiple tests within a single nf-test file, keeping thme in separate files makes it easier to launch individual pipeline tests. Each nf-test file should be named after the scenario it tests in the following format:

- `tests/scenarioA.main.nf.test`
- `tests/scenarioB.main.nf.test`

### Pipeline nf-tests additional guidance

The same guidelines for test profiles, test data and nf-test also apply to pipeline tests. In addition, the following guidelines apply:

- To ensure all output files are caught, the `params.outdir` should be set the the nf-test variable `outputDir`
- The tag `PIPELINE` and the pipeline name should be added to all tests
