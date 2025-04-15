---
title: "nf-test: Writing tests"
subtitle: Guidelines for writing nf-test tests
shortTitle: Writing nf-test tests
weight: 10
---

## Philosophy of nf-tests for nf-core components

- Each component contains a `tests/` folder beside the `main.nf` of the component itself, containing the test files.
- Test files come with a [snapshot](https://code.askimed.com/nf-test/docs/assertions/snapshots/) of component output channels.

## nf-test guidelines for a simple un-chained module

- Some modules MAY require additional parameters added to the test command to successfully run. These can be specified using a params input and an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

If your module requires a `nextflow.config` file to run, create the file to the module's `tests/` directory and add the following code to use parameters defined in the `when` scope of the test.

```bash
touch modules/nf-core/<tool>/<subtool>/tests/nextflow.config
```

```groovy title="nextflow.config"
process {
  withName: 'MODULE' {
    ext.args = params.module_args
  }
}
```

You do not need to modify the contents of this file any further.

Then add the config to the `main.nf.test` file.

```groovy title="main.nf.test"
process "MODULE"
config "./nextflow.config"
```

Lastly supply the params in the when section of the test.

```groovy title="main.nf.test"
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

- When your test data is too big, the tests take too long or require too much resources, you can opt to run your tests in stub mode by adding the following option:

  ```groovy title="main.nf.test"
  options "-stub"
  ```

  :::note
  this can be added at the top of `main.nf.test` to have all tests run in stub mode or this can also be added to a single test
  :::

:::tip
See the [assertions documentation](/docs/contributing/nf-test/assertions) for examples on how to handle different types of test data and scenarios.
:::

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
                            file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)),
                        tuple([ id:'test2', single_end:false ],
                            file(params.modules_testdata_base_path + 'genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz', checkIfExists: true))
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
