# main_nf

#### `SubworkflowLint.main_nf(subworkflow: NFCoreComponent) → Tuple[List[str], List[str]]{:python}`

Lint a `main.nf` subworkflow file

Can also be used to lint local subworkflow files,
in which case failures will be reported as
warnings.

The test checks for the following:

- A subworkflow SHOULD import at least two modules
- All included modules or subworkflows are used and their names are used for versions.yml
- The workflow name is all capital letters
- The subworkflow emits a software version
