# module_tests

#### ModuleLint.module_tests(module)

Lint the tests of a module in `nf-core/modules`

It verifies that the test directory exists
and contains a `main.nf` and a `test.yml`,
and that the module is present in the `pytest_modules.yml`
file.
