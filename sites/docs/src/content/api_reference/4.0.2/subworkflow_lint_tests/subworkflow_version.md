# subworkflow_version

#### `SubworkflowLint.subworkflow_version(subworkflow){:python}`

Verifies that the subworkflow has a version specified in the `modules.json` file

It checks whether the subworkflow has an entry in the `modules.json` file
containing a commit SHA. If that is true, it verifies that there are no
newer version of the subworkflow available.

The following checks are performed:

- `git_sha`: The subworkflow must have a `git_sha` entry in `modules.json`.
- `subworkflow_version`: The subworkflow version must match the latest commit
  in the remote repository. A warning is issued if a newer version is available.
