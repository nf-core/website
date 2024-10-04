# subworkflow_tests

#### `SubworkflowLint.subworkflow_tests(subworkflow: NFCoreComponent){:python}`

Lint the tests of a subworkflow in `nf-core/modules`

It verifies that the test directory exists
and contains a `main.nf.test` a `main.nf.test.snap` and `tags.yml`.

Additionally, hecks that all included components in test `main.nf` are specified in `test.yml`
