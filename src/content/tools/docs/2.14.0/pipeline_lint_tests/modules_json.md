# nextflow_config

#### `PipelineLint.modules_json(){:python}`

Make sure all modules described in the `modules.json` file are actually installed

Every module installed from `nf-core/modules` must have an entry in the `modules.json` file
with an associated version commit hash.

- Failure: If module entries are found in `modules.json` for modules that are not installed
