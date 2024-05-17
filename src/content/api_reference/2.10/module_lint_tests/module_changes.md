# module_changes

#### `ModuleLint.module_changes(module){:python}`

Checks whether installed nf-core modules have changed compared to the
original repository

Downloads the `main.nf` and `meta.yml` files for every module
and compares them to the local copies

If the module has a commit SHA entry in the `modules.json`, the file content is
compared against the files in the remote at this SHA.

Only runs when linting a pipeline, not the modules repository
