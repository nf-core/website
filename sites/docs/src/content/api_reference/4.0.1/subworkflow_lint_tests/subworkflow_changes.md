# subworkflow_changes

#### `SubworkflowLint.subworkflow_changes(subworkflow){:python}`

Checks whether installed nf-core subworkflow have changed compared to the
original repository

Downloads the `main.nf` and `meta.yml` files for every subworkflow
and compares them to the local copies

If the subworkflow has a commit SHA entry in the `modules.json`, the file content is
compared against the files in the remote at this SHA.

Only runs when linting a pipeline, not the modules repository

The following checks are performed:

- `subworkflow_patch`: If the subworkflow is patched, the patch must apply
  cleanly in reverse against the remote version.
- `check_local_copy`: Each subworkflow file must be identical to the
  corresponding file in the remote repository at the pinned commit SHA.
