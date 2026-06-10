# module_tests

#### `ModuleLint.module_tests(module: NFCoreComponent, allow_missing: bool = False){:python}`

Lint the tests of a module in `nf-core/modules`

Checks the `tests/` directory and `main.nf.test` file for correctness,
validates snapshot content, and verifies that nf-test tags follow guidelines.

The following checks are performed:

- `test_dir_exists`: The nf-test directory `tests/` must exist.
- `test_main_nf_exists`: The file `tests/main.nf.test` must exist.
- `test_snapshot_exists`: If `snapshot()` is called in `main.nf.test`,
  the snapshot file `tests/main.nf.test.snap` must exist and be valid JSON.
- `test_snap_md5sum`: The snapshot must not contain md5sums for empty files
  (`d41d8cd98f00b204e9800998ecf8427e`) or empty compressed files
  (`7029066c27ac6f5ef18d660d5741979a`), unless the test name contains `stub`.
- `test_snap_versions`: The snapshot must contain a `versions` key.
- `test_main_tags`: The `main.nf.test` file must declare the required tags:
  `modules`, `modules_<org>`, the full component name, and (for `tool/subtool`
  modules) the tool name alone. Tags for any chained components included via
  `include` statements must also be present.
- `test_old_test_dir`: The legacy pytest directory
  `tests/modules/<component_name>/` must not exist.
