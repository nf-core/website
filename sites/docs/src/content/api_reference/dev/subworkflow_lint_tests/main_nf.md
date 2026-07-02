# main\_nf

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

## main\_nf\_exists

The `main.nf` file must exist.

## main\_nf\_script\_outputs

The workflow must have an `emit:` block.

## main\_nf\_version\_emitted

The subworkflow should emit a software version
channel. A warning is issued if no `versions` output is found (can be
ignored if the subworkflow uses topic channels).
