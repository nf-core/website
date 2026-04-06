# subworkflow_tests

#### `SubworkflowLint.subworkflow_tests(subworkflow: NFCoreComponent, allow_missing: bool = False){:python}`

Lint the tests of a subworkflow in `nf-core/modules`

It verifies that the test directory exists
and contains a `main.nf.test` and a `main.nf.test.snap`

Additionally, checks that all included components in test `main.nf` are specified in `test.yml`
