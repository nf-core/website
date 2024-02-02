# subworkflow_tests

#### SubworkflowLint.subworkflow_tests(subworkflow)

Lint the tests of a subworkflow in `nf-core/modules`

It verifies that the test directory exists
and contains a `main.nf` and a `test.yml`,
and that the subworkflow is present in the `pytest_modules.yml`
file.

Additionally, hecks that all included components in test `main.nf` are specified in `test.yml`
