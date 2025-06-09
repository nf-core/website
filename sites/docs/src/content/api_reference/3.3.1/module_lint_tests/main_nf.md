# main_nf

#### `ModuleLint.main_nf(module: NFCoreComponent, fix_version: bool, registry: str, progress_bar: Progress) → Tuple[List[str], List[str]]{:python}`

Lint a `main.nf` module file

Can also be used to lint local module files,
in which case failures will be reported as
warnings.

The test checks for the following:

- Software versions and containers are valid
- The module has a process label and it is among
  the standard ones.
- If a `meta` map is defined as one of the modules
  inputs it should be defined as one of the outputs,
  and be correctly configured in the `saveAs` function.
- The module script section should contain definitions
  of `software` and `prefix`
