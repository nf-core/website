# subworkflow_tests

#### `SubworkflowLint.subworkflow_tests(subworkflow: NFCoreComponent, allow_missing: bool = False){:python}`

Lint the tests of a subworkflow in `nf-core/modules`

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
- `test_snap_versions`: The snapshot should contain a `versions` key.
  A warning (not a failure) is issued if it is absent, since subworkflows that
  use topic channels may not emit versions directly.
- `test_main_tags`: The `main.nf.test` file must declare the required tags:
  `subworkflows`, `subworkflows/<component_name>`, `subworkflows_<org>`,
  all components included in the subworkflow’s `main.nf`, and any chained
  components referenced via `include` statements in the test file.
- `test_old_test_dir`: The legacy pytest directory
  `tests/subworkflows/<component_name>/` must not exist.
