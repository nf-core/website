# subworkflow_version

#### `SubworkflowLint.subworkflow_version(subworkflow){:python}`

Verifies that the subworkflow has a version specified in the `modules.json` file

It checks whether the subworkflow has an entry in the `modules.json` file
containing a commit SHA. If that is true, it verifies that there are no
newer version of the subworkflow available.
