# module\_tests

#### `ModuleLint.module_tests(module: NFCoreComponent, allow_missing: bool = False){:python}`

Lint the tests of a module in `nf-core/modules`

Checks the `tests/` directory and `main.nf.test` file for correctness,
validates snapshot content, and verifies that nf-test tags follow guidelines.

The following checks are performed:

## test\_dir\_exists

The nf-test directory `tests/` must exist.

## test\_main\_nf\_exists

The file `tests/main.nf.test` must exist.

## test\_snapshot\_exists

If `snapshot()` is called in `main.nf.test`,
the snapshot file `tests/main.nf.test.snap` must exist and be valid JSON.

## test\_snap\_md5sum

The snapshot must not contain md5sums for empty files
(`d41d8cd98f00b204e9800998ecf8427e`) or empty compressed files
(`7029066c27ac6f5ef18d660d5741979a`), unless the test name contains `stub`.

## test\_snap\_versions

The snapshot must contain a `versions` key.

## test\_main\_tags

The `main.nf.test` file must declare the required tags:
`modules`, `modules_<org>`, the full component name, and (for `tool/subtool`
modules) the tool name alone. Tags for any chained components included via
`include` statements must also be present.

## test\_old\_test\_dir

The legacy pytest directory
`tests/modules/<component_name>/` must not exist.

## test\_stub\_gzip\_syntax

Stub blocks with gzip output files must use the proper syntax. Stub files
ending in `.gz` must use `echo "" | gzip > file.gz`. Simply touching or
creating empty `.gz` files will break nf-test’s gzip parser.
