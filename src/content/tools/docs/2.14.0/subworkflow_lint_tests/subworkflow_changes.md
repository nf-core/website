# subworkflow_changes

#### `SubworkflowLint.subworkflow_changes(subworkflow){:python}`

Checks whether installed nf-core subworkflow have changed compared to the
original repository

Downloads the `main.nf` and `meta.yml` files for every subworkflow
and compares them to the local copies

If the subworkflow has a commit SHA entry in the `modules.json`, the file content is
compared against the files in the remote at this SHA.

Only runs when linting a pipeline, not the modules repository
