# main_nf

#### `SubworkflowLint.main_nf(subworkflow: NFCoreComponent) → tuple[list[str], list[str]]{:python}`

Lint a `main.nf` subworkflow file

Can also be used to lint local subworkflow files,
in which case failures will be reported as
warnings.

The test checks for the following:

- A subworkflow SHOULD import at least two modules
- All included modules or subworkflows are used and their names are used for versions.yml
- The workflow name is all capital letters
- The subworkflow emits a software version

The following checks are performed:

- `main_nf_exists`: The `main.nf` file must exist.
- `main_nf_script_outputs`: The workflow must have an `emit:` block.
- `main_nf_version_emitted`: The subworkflow should emit a software version
  channel. A warning is issued if no `versions` output is found (can be
  ignored if the subworkflow uses topic channels).
